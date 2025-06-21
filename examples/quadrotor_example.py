import casadi as ca
import numpy as np

# Physical parameters
m = 1.0
g = 9.81
l = 0.25
kF = 1e-5
kM = 1e-6
Ix, Iy, Iz = 0.01, 0.01, 0.02  # Diagonal inertia terms

# rotor limits
omega_min = 0.0
omega_max = 600.0
# orientation limits
phi_min = -ca.pi / 4
phi_max = ca.pi / 4
theta_min = -ca.pi / 4
theta_max = ca.pi / 4

# State variables
x, y, z = ca.MX.sym('x'), ca.MX.sym('y'), ca.MX.sym('z')
vx, vy, vz = ca.MX.sym('vx'), ca.MX.sym('vy'), ca.MX.sym('vz')
phi, theta, psi = ca.MX.sym('phi'), ca.MX.sym('theta'), ca.MX.sym('psi')
p, q, r = ca.MX.sym('p'), ca.MX.sym('q'), ca.MX.sym('r')

state = ca.vertcat(x, y, z, vx, vy, vz, phi, theta, psi, p, q, r)

# Control inputs
w1, w2, w3, w4 = ca.MX.sym('w1'), ca.MX.sym('w2'), ca.MX.sym('w3'), ca.MX.sym('w4')
omega = ca.vertcat(w1, w2, w3, w4)

# Rotor thrusts
thrusts = kF * omega**2
T = ca.sum1(thrusts)

# Torques
tau_phi = l * kF * (w2**2 - w4**2)
tau_theta = l * kF * (w3**2 - w1**2)
tau_psi = kM * (w1**2 - w2**2 + w3**2 - w4**2)

Rx = ca.vertcat(
    ca.horzcat(1, 0, 0),
    ca.horzcat(0, ca.cos(phi), -ca.sin(phi)),
    ca.horzcat(0, ca.sin(phi), ca.cos(phi))
)

# Rotation matrix around Y-axis (pitch)
Ry = ca.vertcat(
    ca.horzcat(ca.cos(theta), 0, ca.sin(theta)),
    ca.horzcat(0, 1, 0),
    ca.horzcat(-ca.sin(theta), 0, ca.cos(theta))
)

# Rotation matrix around Z-axis (yaw)
Rz = ca.vertcat(
    ca.horzcat(ca.cos(psi), -ca.sin(psi), 0),
    ca.horzcat(ca.sin(psi), ca.cos(psi), 0),
    ca.horzcat(0, 0, 1)
)

# Combined rotation matrix: R = Rz * Ry * Rx
R = Rz @ Ry @ Rx
#  psi thea phi
#  r   q    p

# Translational acceleration
acc = (1/m) * (R.T @ ca.vertcat(0, 0, T)) - ca.vertcat(0, 0, g)

# Euler angle rates
# Transformation matrix from body angular velocities to Euler angle rates
# This matrix is valid for a ZYX (Yaw-Pitch-Roll) Euler angle sequence.
# It relates (p, q, r) to (phi_dot, theta_dot, psi_dot)
euler_angle_transform_matrix = ca.vertcat(
    ca.horzcat(1, ca.sin(phi) * ca.tan(theta), ca.cos(phi) * ca.tan(theta)),
    ca.horzcat(0, ca.cos(phi), -ca.sin(phi)),
    ca.horzcat(0, ca.sin(phi) / ca.cos(theta), ca.cos(phi) / ca.cos(theta))
)

# Vector of body angular velocities
body_angular_velocities = ca.vertcat(p, q, r)

# Calculate the Euler angle derivatives
# euler_dot = [phi_dot, theta_dot, psi_dot]
euler_dot = euler_angle_transform_matrix @ body_angular_velocities

# Angular velocity derivatives (simplified)
p_dot = (1/Ix) * (tau_phi + (Iy - Iz)*q*r)
q_dot = (1/Iy) * (tau_theta + (Iz - Ix)*p*r)
r_dot = (1/Iz) * (tau_psi + (Ix - Iy)*p*q)

# Concatenate derivative
state_dot = ca.vertcat(
    vx, vy, vz,
    acc,
    euler_dot,
    p_dot, q_dot, r_dot
)

# CasADi function
f_dyn = ca.Function('f', [state, omega], [state_dot]).expand()

## constraints

# path constraints

f_control = ca.Function('f_control', [state, omega], [ca.vertcat(omega, phi, theta)]).expand()
lb_control = ca.vertcat(omega_min*ca.DM.ones(4), phi_min, theta_min)
ub_control = ca.vertcat(omega_max*ca.DM.ones(4), phi_max, theta_max)

# t0 #

# Initial state (hover at origin)
x0 = ca.DM([
    0, 0, 0,      # position x, y, z
    0, 0, 0,      # linear velocities
    0, 0, 0,      # roll, pitch, yaw
    0, 0, 0       # angular velocities
])
f_t0 = ca.Function('f_t0', [state], [state - x0]).expand()
lb_t0 = ca.DM.zeros(12)  # Equality constraints for initial state
ub_t0 = ca.DM.zeros(12)  # Equality constraints for initial state 

# tN #

# Final state (target position, hover)
Xn = ca.DM([
    5, 5, 2,      # position
    0, 0, 0,      # linear velocities
    0, 0, 0,      # roll, pitch, yaw
    0, 0, 0       # angular velocities
])

f_tN = ca.Function('f_tN', [state], [state - Xn]).expand()
lb_tN = ca.DM.zeros(12)  # Equality constraints for final state
ub_tN = ca.DM.zeros(12)  # Equality constraints for final state

# cost function
f_cost = ca.Function('f_cost', [omega], [ca.sumsqr(omega/omega_max)]).expand()

# set up the optimization problem

K = 100 
T = 5.0  # Total time

nx = [12 for _ in range(K)]
nu = [4 for _ in range(K-1)] + [0]

f_intg = ca.Function('f_dyn_discr', [state, omega], [state + T/K*f_dyn(state, omega)]).expand()

def discrete_dynamics(uk, xk, k):
    return f_intg(xk, uk)

def cost(uk, xk, k):
    if k < K-1:
        return f_cost(uk)
    return 0.0 # no terminal cost

# expects a list of tuples in the form (lb, f, ub)
def path_constraints(uk, xk ,k):
    cc = []
    # equality constraints
    if k == 0:
        cc.append((lb_t0, f_t0(xk), ub_t0))
        cc.append((lb_control, f_control(xk, uk), ub_control))
    elif k == K-1:
        cc.append((lb_tN, f_tN(xk), ub_tN))
    else:
        cc.append((lb_control, f_control(xk, uk), ub_control))
    return cc

def initial_guess(k):
    if k == K-1:
        u_init = ca.DM(1)
    else:
        u_init = np.sqrt(g * m / (4 * kF)) * np.ones(4) # hover rotor speed
    alpha = k / (K - 1)
    x_init = (1 - alpha) * x0 + alpha * Xn + np.ones(12) * 0.1  # interpolate between start and end state with a small offset to avoid singularities
    return u_init, x_init

opti = ca.Opti()
x = []
u = []
for k in range(K):
    x.append(opti.variable(nx[k]))
    u.append(opti.variable(nu[k]))
# add constraints
ng = []
for k in range(K):
    if (k < K-1):
        opti.subject_to(x[k+1] == discrete_dynamics(u[k], x[k], k))
    path_constr = path_constraints(u[k], x[k], k)
    ng.append(0)
    for constr in path_constr:
        ng[-1] += constr[1].nnz()
        opti.subject_to((constr[0]<=constr[1])<=constr[2])
# set the objective
J = 0
for k in range(K):
    J += cost(u[k], x[k], k)
opti.minimize(J)

# set initial guesses
for k in range(K):
    u_init, x_init = initial_guess(k)
    opti.set_initial(u[k], u_init)
    opti.set_initial(x[k], x_init)


# solve with Ipopt
opti.solver('ipopt', {'ipopt.print_info_string': 'yes'})
opti.solve()

# solve with Fatrop
opti.solver('fatrop', {'structure_detection':'manual', 'nx': nx, 'nu':nu, 'ng':ng, 'N':K-1, "expand": True, "fatrop.mu_init":1e-1, "jit":False, "fatrop.print_level":10, "jit_options": {"flags": "-O3", "verbose": True}})
res = opti.solve()

X_sol = res.value(ca.horzcat(*x))
print("Solution found:", res.stats()['success'])

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Assuming X is your [12, N] state array
x = X_sol[0, :]
y = X_sol[1, :]
z = X_sol[2, :]

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Trajectory
ax.plot(x, y, z, label='Drone Trajectory', color='blue', linewidth=2)
ax.scatter(x[0], y[0], z[0], color='green', label='Start', s=50)
ax.scatter(x[-1], y[-1], z[-1], color='red', label='End', s=50)

# Equal axis scaling workaround
max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
mid_x = (x.max()+x.min()) * 0.5
mid_y = (y.max()+y.min()) * 0.5
mid_z = (z.max()+z.min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Labels and formatting
ax.set_title('3D Drone Flight Trajectory')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.legend()
ax.grid(True)
ax.view_init(elev=30, azim=135)

plt.tight_layout()
plt.show()