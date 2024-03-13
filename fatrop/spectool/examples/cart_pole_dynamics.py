import casadi as cs
import numpy as np


m_cart = 1.0
m_pole = 0.1
g = 9.81
l = 0.5
I = m_pole * l**2 / 3

def ode(x_cart, dx_cart, theta_pole, dtheta_pole, F):
    #  helper functions
    def time_der(expr):
        return [cs.jacobian(expr, q), cs.jacobian(expr, dq)] #  @ [dq, ddq]
    def solve_2x2(A, b):
        def det_2x2(A):
            return A[0,0]*A[1,1] - A[0,1]*A[1,0]
        detA = det_2x2(A)
        x0  = det_2x2(cs.horzcat(b, A[:,1]))/detA
        x1  = det_2x2(cs.horzcat(A[:,0], b))/detA
        return cs.vertcat(x0, x1)

    q = cs.vertcat(x_cart, theta_pole)
    dq = cs.vertcat(dx_cart, dtheta_pole)
    
    # using Lagrangian dynamics
    T = 0.5 * m_cart * dx_cart**2 + 0.5 * m_pole * (dx_cart + l/2 * dtheta_pole)**2
    V = m_pole * g * l/2 * cs.cos(theta_pole)
    W = F * x_cart
    L = T - V + W


    dt_L_dq = time_der(cs.gradient(L, dq))
    A = dt_L_dq[1]
    b = cs.gradient(L, q) - dt_L_dq[0]@dq
    ddq = solve_2x2(A, b)
    return dx_cart, ddq[0], dtheta_pole, ddq[1]

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