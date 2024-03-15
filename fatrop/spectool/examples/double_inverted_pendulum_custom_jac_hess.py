import casadi as cs 
import numpy as np
import fatropy.spectool as sp
import double_inverted_pendulum_dynamics
import robot2d as r2d

# control grid parameters
N = 50 # number of time steps
T = 2.0 # time horizon

# Instantiate an Ocp 
ocp = sp.Ocp()

# add a stage to the ocp (N = 50 time steps)
stage = ocp.new_stage(N)

# define the state and control variables
theta_pole = ocp.state(2) # angle of the pole, measured from the vertical
dtheta_pole = ocp.state(2) # angular velocity of the pole
F = ocp.control() # horizontal force applied to the cart
x = [theta_pole, dtheta_pole]
ux = [F] + x

# create your problem parameters -> value at time step 0 for each state
theta_pole_0 = ocp.parameter(2)
dtheta_pole_0 = ocp.parameter(2)

# Define the system dynamics
x_dot, mechanism, jac, hess = double_inverted_pendulum_dynamics.ode(*x, F)

# # take a runge kutta 4 integrator
# intg = sp.IntegratorRk4(list(zip(x, x_dot)), T/N)
# for x, xn, in zip(x, intg(x)):
#     stage.set_next(x, xn) 
# or a simple euler integrator with a custom Jacbian
stage.set_next(theta_pole, theta_pole + T/N*x_dot[0])
jacobian_struct = sp.Jacobian(cs.jacobian(dtheta_pole, cs.vcat(ux))+  T/N*jac, cs.vcat(ux))
print((T/N*hess[1]).shape)
print(cs.gradient(dtheta_pole.T@ hess[2], cs.vcat(ux)).shape)
hessian_struct = sp.Hessian(T/N*hess[0], cs.gradient(dtheta_pole.T@ hess[2], cs.vcat(ux)) + T/N*hess[1], cs.vcat(ux), hess[2])
stage.set_next(dtheta_pole, dtheta_pole + T/N*x_dot[1],
               jacobian =  jacobian_struct, 
               hessian  =  hessian_struct)
# stage.set_next(dtheta_pole, dtheta_pole + T/N*x_dot[1])

# Define the cost function
h_pole = double_inverted_pendulum_dynamics.l/2*(cs.cos(theta_pole[0]) + cs.cos(theta_pole[1]))
stage.add_objective(-100*h_pole, sp.t0, sp.mid, sp.tf) # we want the pole to be upright
stage.add_objective(cs.sumsqr(dtheta_pole), sp.t0, sp.mid, sp.tf) # we want the pole to be upright

# Define the constraints
stage.at_t0().subject_to(theta_pole == theta_pole_0)
stage.at_t0().subject_to(dtheta_pole == dtheta_pole_0)

# pick a solver
ocp.solver("fatrop", {'post_expand':True, 'jit':True}, {}) # fatrop/ipopt

# put the solver inside a function
func = ocp.to_function("double_pendulum", [theta_pole_0, dtheta_pole_0], [ocp.at_t0(F), ocp.sample(theta_pole)[1]])

# solve the ocp
u0_sol, theta_pol_sol = func(np.array([0.1, 0.0]), np.array([0.0, 0.0]))
r2d.animate(((mechanism, theta_pole, theta_pol_sol),))