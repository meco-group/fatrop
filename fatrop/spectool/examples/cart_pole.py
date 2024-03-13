import casadi as cs 
import numpy as np
import fatropy.spectool as sp
import cart_pole_dynamics

# control grid parameters
N = 50 # number of time steps
T = 2.0 # time horizon

# Instantiate an Ocp 
ocp = sp.Ocp()

# add a stage to the ocp (N = 50 time steps)
stage = ocp.new_stage(N)

# define the state and control variables
x_cart = ocp.state() # horizontal position of the cart
dx_cart = ocp.state() # horizontal velocity of the cart
theta_pole = ocp.state() # angle of the pole, measured from the vertical
dtheta_pole = ocp.state() # angular velocity of the pole
F = ocp.control() # horizontal force applied to the cart
x = [x_cart, dx_cart, theta_pole, dtheta_pole]

# create your problem parameters -> value at time step 0 for each state
x_cart_0 = ocp.parameter()
dx_cart_0 = ocp.parameter()
theta_pole_0 = ocp.parameter()
dtheta_pole_0 = ocp.parameter()

# Define the system dynamics
x_dot = cart_pole_dynamics.ode(*x, F)

# take a runge kutta 4 integrator
intg = sp.IntegratorRk4(list(zip(x, x_dot)), T/N)
for x, xn, in zip(x, intg(x)):
    stage.set_next(x, xn) 

# Define the cost function
h_pole = cart_pole_dynamics.l/2*cs.cos(theta_pole)
stage.add_objective(-h_pole, sp.t0, sp.mid, sp.tf)
stage.add_objective(dx_cart**2, sp.t0, sp.mid, sp.tf)

# Define the constraints
stage.at_t0().subject_to(x_cart == x_cart_0)
stage.at_t0().subject_to(dx_cart == dx_cart_0)
stage.at_t0().subject_to(theta_pole == theta_pole_0)
stage.at_t0().subject_to(dtheta_pole == dtheta_pole_0)
stage.subject_to((-1. <= x_cart) <= 1., sp.t0, sp.mid, sp.tf)

# pick a solver
ocp.solver("fatrop", {'post_expand':True, 'jit':True}, {}) # fatrop/ipopt

# put the solver inside a function
func = ocp.to_function("cart_pole", [x_cart_0, dx_cart_0, theta_pole_0, dtheta_pole_0], [ocp.at_t0(F), ocp.sample(x_cart)[1], ocp.sample(theta_pole)[1]])

# solve the ocp
u0_sol, x_cart_sol, theta_pol_sol = func(0.0, 0.0, 1.0, 1.0)
print(x_cart_sol)
cart_pole_dynamics.animate(x_cart_sol, theta_pol_sol)