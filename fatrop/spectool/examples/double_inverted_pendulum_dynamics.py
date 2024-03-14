import casadi as cs
import numpy as np
import lagrangian_dynamics
import robot2d as r2d


m_cart = 1.0
m_pole = 1.0
g = 9.81
l = 0.5
I = m_pole * l**2 / 3

def ode(theta_pole, dtheta_pole, F):
    q = cs.vertcat(theta_pole)
    dq = cs.vertcat(dtheta_pole)
    revol1 = r2d.revolute(r2d.transform2d(0, 0, np.pi/2), theta_pole[0], dtheta_pole[0])
    pole1 = r2d.link(revol1, m_pole, l)
    revol2 = r2d.revolute(pole1.right, theta_pole[1], dtheta_pole[1])
    pole2 = r2d.link(revol2, m_pole, l)
    mechanism = r2d.mechanism([pole1, pole2])
    W = F * revol1.theta 
    ddq, jac, hess = mechanism.get_dynamics(W, q, dq, F)
    return [dtheta_pole, ddq], mechanism, jac, hess

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle

def animate(theta_pole_sol):
    theta_pole_sol = np.array(theta_pole_sol) 
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
    ax.set_aspect('equal')
    ax.grid()

    pole, = ax.plot([], [], 'r', lw=2)
    pole2, = ax.plot([], [], 'r', lw=2)

    def init():
        pole.set_data([], [])
        pole2.set_data([], [])
        return pole, pole2

    def animate(i):
        theta = theta_pole_sol[0, i]
        pole.set_data([0, l*cs.sin(theta)], [0, l*cs.cos(theta)])
        pole2.set_data([l*cs.sin(theta), l*cs.sin(theta) + l*cs.sin(theta + theta_pole_sol[1, i])], [l*cs.cos(theta), l*cs.cos(theta) + l*cs.cos(theta + theta_pole_sol[1, i])])
        return pole, pole2

    ani = animation.FuncAnimation(fig, animate, frames=theta_pole_sol.shape[1], init_func=init, blit=True)
    plt.show()