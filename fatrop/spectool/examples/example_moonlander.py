import fatropy.spectool as sp
import numpy as np
import casadi as cs
from moonlander.helpers import transf


# Define the problem parameters (here hard-coded, you can use ocp.register_parameter if you want to be able to change them, cf. below for an example ("target"))
m = 1.0
g = 9.81
I = 0.1
D = 1.0
max_thrust = 2*g*m

# Create an OCP object
ocp = sp.Ocp()
# create a parameter, the name "target" can be used to retrieve the parameter value from the generated interface (see interface_usage.(cpp/py))
target = ocp.parameter(2)

# define the states
p = ocp.state(2)
dp = ocp.state(2)
theta = ocp.state()
dtheta = ocp.state()

# Define the controls
F1 = ocp.control()
F2 = ocp.control()

# define the hybrids (T is the total time here, it will be used in the integrator)
T = ocp.hybrid()


# define the stage with 50 control intervals
stage = ocp.new_stage(50)

# Compute the tranformation of the local frame attached to the moonlander
F_r = transf(theta, p)
F_tot = (F_r @ cs.vertcat(0, F1 + F2, 0)) [:2]

# set up the integrator
intg = sp.IntegratorRk4([
(p, dp),
(theta, dtheta),
(dp, 1/m * F_tot + cs.vertcat(0, -g)),
(dtheta, 1/I * D/2 * (F2 - F1))
    ], T/50.)

# set up the discretized dynamics
stage.set_next(p, intg(p))
stage.set_next(theta, intg(theta))
stage.set_next(dp, intg(dp))
stage.set_next(dtheta, intg(dtheta))
stage.set_next(T, T)

# Define the path constraints
stage.subject_to((0<F1)<max_thrust, sp.t0, sp.mid)
stage.subject_to((0<F2)<max_thrust, sp.t0, sp.mid)

# Define the initial constraints
ocp.at_t0().subject_to(p == [0, 0])
ocp.at_t0().subject_to(dp == [0, 0])
ocp.at_t0().subject_to(theta == 0)
ocp.at_t0().subject_to(dtheta == 0)

# Define the terminal constraints
ocp.at_tf().subject_to(p == target)
ocp.at_tf().subject_to(dp == [0, 0])

# Define the cost function (minimal total time)
ocp.at_tf().add_objective(T)
ocp.at_t0().subject_to(T > 0)



# Provide an initial guess, here we choose a constant over the whole control horizon. It is also possible to provide a nx x N array here.
ocp.set_initial(F1, cs.MX(5.))
ocp.set_initial(F2, cs.MX(5.))
ocp.set_initial(T, cs.MX(10.))

# Retrieve the results
k_F_r_sample, F_r_sample = ocp.sample(F_r)
k_F_tot_sample, F_tot_sample = ocp.sample(F_tot)
k_F1_sample, F1_sample = ocp.sample(F1)
k_F2_sample, F2_sample = ocp.sample(F2)
k_p_sample, p_sample = ocp.sample(p)

### More efficient usage of the solver: make use of to_function, to reduce overhead of the rockit layer
# the arguments should consist of all defined problem parameters
# state and control variables are treated as initial guess for the solver
# stage and control variables that are not defined as arguments use the default initialization
# fatrop_func = ocp.to_function([target], [p_sample, F1_sample, F2_sample, ocp.sample(T)[1], cs.horzcat(0, ocp.sample(T/50.)[1])[:-1]])
ocp.solver('fatrop', {"jit":True}, {"mu_init":1e-1})
fatrop_func = ocp.to_function("moonlander_ocp", [target], [p_sample, F1_sample, F2_sample, ocp.sample(T)[1], cs.horzcat(0, ocp.sample(T/50.)[1])[:-1], ocp.all_variables()])
p_sol, F1_sol, F2_sol, T_sol, timesteps, all_vars = fatrop_func(np.array([6., 5.]))
t_F1_sol = timesteps[k_F1_sample]
print("t_F1_sol", t_F1_sol)
print("F1_sol", F1_sol)
ocp.solver('fatrop', {"jit":True}, {"mu_init":1e-4})
fatrop_func_warmstart = ocp.to_function("moonlander_ocp_warmstart", [target, ocp.all_variables()], [p_sample])
fatrop_func_warmstart(np.array([10., 5.]), all_vars)

# plot the trajectory
import matplotlib.pyplot as plt
plt.plot(p_sol[0, :], p_sol[1,:], 'o-')
plt.show()
