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
        inertia = mass * length**2 / 12
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
    def get_dynamics(self, W, q, qd):
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
        return ld.get_ddq(q, qd, T, V, W)
    
    def animate(self, q_sym, q_vals):
        # collect all symbols for left and right position of each link
        no_links = len(self.links)
        left = []
        right = []
        for link in self.links:
            left.append(link.left.origin())
            right.append(link.right.origin())
        # make a casadi function that evaliates the left and right position of each link
        f = cs.Function('f', [q_sym], [cs.horzcat(*left), cs.horzcat(*right)])

        fig = plt.figure()
        ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
        ax.set_aspect('equal')
        ax.grid()

        axs_link = [ax.plot([], [], 'r', lw=2)[0] for _ in range(no_links)]

        def init():
            for ax_link in axs_link:
                ax_link.set_data([], []) 
            return axs_link 

        def animate(i):
            lefts, rights = f(q_vals[:, i])
            lefts = np.array(lefts)
            rights = np.array(rights)
            for j in range(no_links):
                axs_link[j].set_data([lefts[0, j], rights[0, j]], [lefts[1, j], rights[1, j]])
            return axs_link 

        ani = animation.FuncAnimation(fig, animate, frames=q_vals.shape[1], init_func=init, blit=True)
        plt.show()
            