import casadi as cs
import numpy as np

K = 100 
m = 1.
dt = 0.05

nx = [4 for _ in range(K)]
nu = [2 for _ in range(K-1)] + [0]

def discrete_dynamics(uk, xk, k):
    dp = xk[2:]
    dv = uk / m 
    dv[0] += 5.*uk[1]*uk[1]
    return xk + dt*cs.vertcat(dp, dv)

def cost(uk, xk, k):
    return cs.sumsqr(uk)

def path_constraints(uk, xk ,k):
    cc = []
    # equality constraints
    if k == 0:
        cc.append((0., xk[0] - 0.0, 0))
        cc.append((0., xk[1] - 0.0, 0))
        cc.append((0., xk[2] - 0.0, 0))
        cc.append((0., xk[3] - 0.0, 0))
    elif k == K-1:
        cc.append((0., xk[0] - 10., 0.))
        cc.append((0., xk[1] - 20., 0.))
        cc.append((0., xk[2] - 30., 0.))
        cc.append((0., xk[3] - 40., 0.))

    # inequality constraints
    if k < K-1:
        cc.append((-50. , uk[0], 50.))
        cc.append((-100., uk[1], 100.))
    # no collision
    obstacle_pos = cs.vertcat(5., 10.)
    obstacle_radius = 1. 
    if k < K-1:
        cc.append((0, cs.sumsqr(xk[:2] - obstacle_pos) - obstacle_radius*obstacle_radius, np.inf))
    return cc



opti = cs.Opti()
x = []
u = []
for k in range(K):
    x.append(opti.variable(nx[k]))
    u.append(opti.variable(nu[k]))
# add constraints
for k in range(K):
    if (k < K-1):
        opti.subject_to(x[k+1] == discrete_dynamics(u[k], x[k], k))
    path_constr = path_constraints(u[k], x[k], k)
    for constr in path_constr:
        opti.subject_to((constr[0]<= constr[1])<= constr[2])
# set the objective
J = 0
for k in range(K):
    J += cost(u[k], x[k], k)
opti.minimize(J)

# # opti.to_function("opti_func", [], [opti.x]).generate('casadi_generated.c', {"with_header": True})
# opti.solver('ipopt', {'ipopt.print_info_string': 'yes', 'ipopt.kappa_d':0.})
# # opti.solver('fatrop', {'structure_detection': 'auto', 'fatrop.mu_init': 1e-1})
# opti.solve()
# opti.solver('fatrop', {'structure_detection': 'auto', 'fatrop.mu_init': 1e-1})
# opti.solve()



def l1_pen_problem():
    rho = 1000.
    zeta = 1e-2
    opti = cs.Opti()
    x = []
    u = []
    for k in range(K):
        x.append(opti.variable(nx[k]))
        u.append(opti.variable(nu[k]))
    # add constraints
    p = []
    n = []
    for k in range(K):
        if (k < K-1):
            opti.subject_to(x[k+1] == discrete_dynamics(u[k], x[k], k))
        path_constr = path_constraints(u[k], x[k], k)
        for constr in path_constr:
            pi = opti.variable()
            ni = opti.variable()
            p.append(pi)
            n.append(ni)
            opti.subject_to((constr[0]<= constr[1]+pi-ni)<= constr[2])
            opti.subject_to(pi >= 0)
            opti.subject_to(ni >= 0)
    # set the objective
    J = rho * cs.sum1(cs.vertcat(*p,*n)) + cs.sumsqr(zeta*cs.vertcat(*x, *u))
    opti.minimize(J)
    opti.solver('ipopt')
    opti.solve()
    return opti

opti = l1_pen_problem()

