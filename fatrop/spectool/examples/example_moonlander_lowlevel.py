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

#########
# Create an OCP object
#########
ocp = sp.Ocp()
# create a parameter, the name "target" can be used to retrieve the parameter value from the generated interface (see interface_usage.(cpp/py))
target = ocp.parameter(2)

#########
# define the states
#########
p = ocp.state(2)
dp = ocp.state(2)
theta = ocp.state()
dtheta = ocp.state()

#########
# Define the controls
#########
F1 = ocp.control()
F2 = ocp.control()

#########
# define the T state (T is the total time here, it will be used in the integrator)
#########
T = ocp.state()


#########
# start from a prototype microstage that contains all the states and parameters
#########
ustage_proto = sp.uStage()
ustage_proto.register_state([p, dp, theta, dtheta, T])
ustage_proto.register_global_parameter([target])
#########
# define a proto microstage that contains the dynamics
#########
ustage_proto_dyn = ustage_proto.clone()
ustage_proto_dyn.register_control([F1, F2])
#########
# define the dynamics and the constraints related to the dynamics
#########

# compute the local frame attached to the moonlander
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
ustage_proto_dyn.set_next(p, intg(p))
ustage_proto_dyn.set_next(theta, intg(theta))
ustage_proto_dyn.set_next(dp, intg(dp))
ustage_proto_dyn.set_next(dtheta, intg(dtheta))
ustage_proto_dyn.set_next(T, T)

# Define the path constraints
ustage_proto_dyn.subject_to((0<F1)<max_thrust)
ustage_proto_dyn.subject_to((0<F2)<max_thrust)

#########
# Specify the different microstages of the ocp
#########
ustage_t0 = ustage_proto_dyn.clone()
ustage_mid = ustage_proto_dyn.clone()
ustage_tf = ustage_proto.clone() # no dynamics at the terminal stage

# Define the initial constraints
ustage_t0.subject_to(p == [0, 0])
ustage_t0.subject_to(dp == [0, 0])
ustage_t0.subject_to(theta == 0)
ustage_t0.subject_to(dtheta == 0)

# Define the terminal constraints
ustage_tf.subject_to(p == target)
ustage_tf.subject_to(dp == [0, 0])

# Define the cost function (minimal total time)
ustage_tf.add_objective(T)
ustage_t0.subject_to(T > 0)

#########
# build up the ocp from the microstages
#########
ocp.add_ustage(ustage_t0)
for i in range(49):
    ocp.add_ustage(ustage_mid.clone())
ocp.add_ustage(ustage_tf)


#########
# Provide an initial guess, here we choose a constant over the whole control horizon. It is also possible to provide a nx x N array here.
#########
ocp.set_initial(F1, cs.MX(5.))
ocp.set_initial(F2, cs.MX(5.))
ocp.set_initial(T, cs.MX(10.))

# # Retrieve the results
k_F_r_sample, F_r_sample = ocp.sample(F_r)
k_F_tot_sample, F_tot_sample = ocp.sample(F_tot)
k_F1_sample, F1_sample = ocp.sample(F1)
k_F2_sample, F2_sample = ocp.sample(F2)
k_p_sample, p_sample = ocp.sample(p)

### More efficient usage of the solver: make use of to_function, to reduce overhead of the rockit layer
# the arguments should consist of all defined problem parameters
# state and control variables are treated as initial guess for the solver
# stage and control variables that are not defined as arguments use the default initialization
ocp.solver('fatrop', {"jit":True}, {"mu_init":1e-1})
fatrop_func = ocp.to_function("moonlander_ocp", [target], [p_sample, F1_sample, F2_sample, ocp.sample(T)[1], cs.horzcat(0, ocp.sample(T/50.)[1])[:-1], ocp.all_variables()])
p_sol, F1_sol, F2_sol, T_sol, timesteps, all_vars = fatrop_func(np.array([6., 5.]))
t_F1_sol = timesteps[k_F1_sample]
print("t_F1_sol", t_F1_sol)
print("F1_sol", F1_sol)
fatrop_func_warmstart = ocp.to_function("moonlander_ocp_warmstart", [target, ocp.all_variables()], [p_sample])
fatrop_func_warmstart(np.array([10., 5.]), all_vars)

# plot the trajectory
import matplotlib.pyplot as plt
plt.plot(p_sol[0, :], p_sol[1,:], 'o-')
plt.show()
