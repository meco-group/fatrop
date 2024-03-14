import casadi as cs
import numpy as np
import lagrangian_dynamics
import robot2d as r2d


m_cart = 1.0
m_pole = 1.0
g = 9.81
l = 0.5
I = m_pole * l**2 / 3

def ode(x_cart, dx_cart, theta_pole, dtheta_pole, F):
    q = cs.vertcat(x_cart, theta_pole)
    dq = cs.vertcat(dx_cart, dtheta_pole)
    prism = r2d.prismatic(r2d.zero_transform(), x_cart, dx_cart, 0, 0)
    cart = r2d.link(prism, m_cart, 0)
    revol = r2d.revolute(cart.center.transform(r2d.transform2d(0, 0, np.pi/2)), theta_pole, dtheta_pole)
    pole = r2d.link(revol, m_pole, l)
    mechanism = r2d.mechanism([cart, pole])
    W = F * cart.center.x
    ddq, _, _ = mechanism.get_dynamics(W, q, dq, F)
    
    # # # using Lagrangian dynamics
    # T = 0.5 * m_cart * dx_cart**2 + 0.5 * m_pole * (dx_cart + l/2 * dtheta_pole)**2
    # V = m_pole * g * l/2 * cs.cos(theta_pole)
    # L = T - V + W
    return dx_cart, ddq[0], dtheta_pole, ddq[1], mechanism

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle

def animate(x_cart_sol, theta_pole_sol):
    x_cart_sol = np.array(x_cart_sol)
    theta_pole_sol = np.array(theta_pole_sol) 
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
    ax.set_aspect('equal')
    ax.grid()

    pole, = ax.plot([], [], 'r', lw=2)

    def init():
        pole.set_data([], [])
        return pole,

    def animate(i):
        x = x_cart_sol[0, i]
        theta = theta_pole_sol[0, i]
        pole.set_data([x, x + l*cs.sin(theta)], [0, l*cs.cos(theta)])
        return pole,

    ani = animation.FuncAnimation(fig, animate, frames=x_cart_sol.shape[1], init_func=init, blit=True)
    plt.show()