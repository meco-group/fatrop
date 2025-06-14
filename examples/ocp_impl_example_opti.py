import casadi as cs

K = 100 
m = 1.
dt = 0.05

nx = [4 for _ in range(K)]
nu = [2 for _ in range(K-1)] + [0]

def discrete_dynamics(uk, xk, k):
    dp = xk[2:]
    dv = uk / m 
    dv[0] += .5*uk[1]*uk[1]
    return xk + dt*cs.vertcat(dp, dv)

def cost(uk, xk, k):
    return cs.sumsqr(uk)

def path_constraints(uk, xk ,k):
    cc = []
    # equality constraints
    if k == 0:
        cc.append(xk[0] - 0.0 == 0.)
        cc.append(xk[1] - 0.0 == 0.)
        cc.append(xk[2] - 0.0 == 0.)
        cc.append(xk[3] - 0.0 == 0.)
    elif k == K-1:
        cc.append(xk[0] - 1.0 == 0.)
        cc.append(xk[1] - 2.0 == 0.)
        cc.append(xk[2] - 3.0 == 0.)
        cc.append(xk[3] - 4.0 == 0.)

    # inequality constraints
    if k < K-1:
        cc.append(-50. <=(uk[0] <= 50.))
        cc.append(-100.<=(uk[1]<=100.))
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
        opti.subject_to(constr)
# set the objective
J = 0
for k in range(K):
    J += cost(u[k], x[k], k)
opti.minimize(J)

opti.solver('fatrop', {'structure_detection':'manual', 'nx': nx, 'nu':nu, 'ng':ng, 'N':K-1, "expand": True, "fatrop.mu_init":1e-1, "jit":False, "jit_options": {"flags": "-O3", "verbose": True}})
opti.to_function("opti_func", [], [opti.x]).generate('casadi_generated.c', {"with_header": True})
# opti.solver('ipopt', {'ipopt.print_info_string': 'yes', 'ipopt.kappa_d':0.})
opti.solve()