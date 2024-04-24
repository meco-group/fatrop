import casadi as cs
import numpy as np
import lagrangian_dynamics as ld
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class transform2d:
    def __init__(self, x:cs.MX, y:cs.MX, theta:cs.MX, dx:cs.MX = 0., dy:cs.MX =0., dtheta:cs.MX = 0.):
        self.x = x
        self.y = y
        self.theta = theta
        self.dx = dx
        self.dy = dy
        self.dtheta = dtheta
    
    def transformm(self, x, y, theta, dx, dy, dtheta):
        x_ = x*cs.cos(self.theta) - y*cs.sin(self.theta) + self.x
        y_ = x*cs.sin(self.theta) + y*cs.cos(self.theta) + self.y
        theta_ = theta + self.theta
        dx_ = dx*cs.cos(self.theta) - dy*cs.sin(self.theta) - self.dtheta*x*cs.sin(self.theta) - self.dtheta*y*cs.cos(self.theta) + self.dx
        dy_ = dx*cs.sin(self.theta) + dy*cs.cos(self.theta) + self.dtheta*x*cs.cos(self.theta) - self.dtheta*y*cs.sin(self.theta) + self.dy
        dtheta_ = dtheta + self.dtheta
        return transform2d(x_, y_, theta_, dx_, dy_, dtheta_)
    
    def origin(self):
        return cs.vertcat(self.x, self.y)

    def transform(self, other):
        return self.transformm(other.x, other.y, other.theta, other.dx, other.dy, other.dtheta)
    
class zero_transform(transform2d):
    def __init__(self):
        super().__init__(0., 0., 0., 0., 0., 0.)

class rigid_body:
    def __init__(self, mass_center:transform2d, mass:cs.MX, inertia:cs.MX):
        self.mass_center = mass_center
        self.mass = mass
        self.inertia = inertia
    
class joint:
    def __init__(self, pos:transform2d, theta, dtheta, x, dx, y, dy):
        self.theta = theta
        self.dtheta = dtheta
        self.pos = pos
        self.x = x
        self.dx = dx
        self.y = y
        self.dy = dy
    
class prismatic(joint):
    def __init__(self, pos:transform2d, x, dx, y, dy):
        super().__init__(pos, 0, 0, x, dx, y, dy)

class revolute(joint):
    def __init__(self, pos:transform2d, theta, dtheta):
        super().__init__(pos, theta, dtheta, 0, 0, 0, 0)

class world(joint):
    def __init__(self):
        super().__init__(transform2d(0., 0., 0.), 0., 0., 0., 0., 0., 0.)

class link(rigid_body):
    def __init__(self, left_joint:joint, mass:cs.MX, length:cs.MX):
        inertia = mass * np.abs(length)**2 / 12
        self.length = length
        self.left_joint = left_joint
        self.left = self.left_joint.pos.transform(transform2d(self.left_joint.x, self.left_joint.y, self.left_joint.theta, self.left_joint.dx, self.left_joint.dy, self.left_joint.dtheta))
        self.right = self.left.transform(transform2d(self.length, 0, 0))
        self.center = self.left.transform(transform2d(self.length/2, 0, 0))
        super().__init__(self.center, mass, inertia)

# a mechanism is a collection of links
class mechanism:
    def __init__(self, links):
        self.links = links
    def get_dynamics(self, W, q, qd, U):
        T = 0
        V = 0
        for link in self.links:
            mass_center_link = link.mass_center
            y = mass_center_link.y
            dx = mass_center_link.dx
            dy = mass_center_link.dy
            dtheta = mass_center_link.dtheta
            T += 0.5 * link.mass * (dx**2 + dy**2) + 0.5 * link.inertia * dtheta**2
            V += link.mass * 9.81 * y
        return ld.get_ddq_jac(q, qd, U, T, V, W)

    @property
    def mass_total(self):
        return sum([link.mass for link in self.links])

    @property
    def center_of_mass(self):
        return sum([cs.vertcat(link.center.x, link.center.y)*link.mass for link in self.links]) / self.mass_total

    @property
    def kinetic_energy(self):
        T = 0
        for link in self.links:
            mass_center_link = link.center
            dx = mass_center_link.dx
            dy = mass_center_link.dy
            dtheta = mass_center_link.dtheta
            T += 0.5 * link.mass * (dx*dx + dy*dy) + 0.5 * link.inertia * dtheta*dtheta
        return T

    @property
    def potential_energy(self):
        V = 0
        for link in self.links:
            mass_center_link = link.center
            y = mass_center_link.y
            V += link.mass * 9.81 * y
        return V

    @property
    def energy(self):
        return self.kinetic_energy + self.potential_energy
    




    
def animate(mechs_qs_vals):
    mechanisms, q_syms, q_valss = zip(*mechs_qs_vals)
    print(mechanisms)
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
    ax.set_aspect('equal')
    ax.grid()
    axs_links = [None for _ in mechanisms]
    fs = [None for _ in mechanisms]
    offsets = [0 for _ in mechanisms]
    pointers = []

    for mechanism, q_sym, q_vals, i in zip(mechanisms, q_syms, q_valss, range(len(mechanisms))):
        pointers += [i for _ in range(q_vals.shape[1])]
        offsets[i] = 0 if i == 0 else offsets[i-1] + q_vals.shape[1]
        # collect all symbols for left and right position of each link
        no_links = len(mechanism.links)
        left = []
        right = []
        for link in mechanism.links:
            left.append(link.left.origin())
            right.append(link.right.origin())
        # make a casadi function that evaliates the left and right position of each link
        fs[i] = cs.Function('f', [q_sym], [cs.horzcat(*left), cs.horzcat(*right)])
        axs_links[i] = [ax.plot([], [], 'r', lw=2)[0] for _ in range(no_links)]

    def init():
        for ax_linksi in axs_links:
            for ax_link in ax_linksi:
                ax_link.set_data([], []) 
        return [ax for axlinski in axs_links for ax in axlinski]

    def animate(i):
        # find the right mechanism
        mechi = pointers[i]
        mechanism = mechanisms[mechi]
        lefts, rights = fs[mechi](q_valss[mechi][:, i-offsets[mechi]])
        lefts = np.array(lefts)
        rights = np.array(rights)
        for j in range(len(mechanism.links)):
            axs_links[mechi][j].set_data([lefts[0, j], rights[0, j]], [lefts[1, j], rights[1, j]])
        return axs_links[mechi] 

    ani = animation.FuncAnimation(fig, animate, frames=sum([q_vals.shape[1] for q_vals in q_valss]), init_func=init, blit=True)
    plt.show()
            