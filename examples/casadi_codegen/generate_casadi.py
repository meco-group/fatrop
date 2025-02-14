import casadi as cs
import numpy as np

# parameters
K = 100 
m = 1.
dt = 0.05

# set up the casadi functions that we will re-used at every time step of the ocp

def dyn_fun():
    uk = cs.MX.sym('uk', 2)
    xk = cs.MX.sym('xk', 4)

    dp = xk[2:]
    dv = uk / m 
    dv[0] += 0.5*uk[1]*uk[1]

    return cs.Function("dyn", [uk, xk],  [xk + dt*cs.vertcat(dp, dv)])

# It is always a good idea to put code that is re-used in a casadi Function, which is expanded for efficiency (changing an MX to an SX Function).
fun = dyn_fun().expand()


def path_constraints_fun():
    xk = cs.MX.sym('xk', 4)
    uk = cs.MX.sym('uk', 2)
    cc = []
    cc.append((-50. , uk[0], 50.))
    cc.append((-100., uk[1], 100.))
    # no collision
    lower, con, upper = zip(*cc)
    return [lower, cs.Function("path_constraints", [uk, xk], [cs.vertcat(*con)]).expand(), upper]

pc_fun = path_constraints_fun()
pc_fun[1] = pc_fun[1].expand()
    

# define our optimal control problem

nx = [4 for _ in range(K)]
nu = [2 for _ in range(K-1)] + [0]
    
def discrete_dynamics(uk, xk, k):
    return fun(uk, xk)

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
        cc.append((0., xk[0] - 1., 0.))
        cc.append((0., xk[1] - 2., 0.))
        cc.append((0., xk[2] - 3., 0.))
        cc.append((0., xk[3] - 4., 0.))
    else:
        lower = pc_fun[0]
        con = pc_fun[1]
        upper = pc_fun[2]
        cc.append((lower, con(uk, xk), upper))
    return cc



opti = cs.Opti()
x = []
u = []
for k in range(K):
    x.append(opti.variable(nx[k]))
    u.append(opti.variable(nu[k]))
# add constraints
ng = []
for k in range(K):
    if (k < K-1):
        opti.subject_to(x[k+1] ==discrete_dynamics(u[k], x[k], k))
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

opti.solver('fatrop', {'structure_detection':'manual', 'nx': nx, 'nu':nu, 'ng':ng, 'N':K-1, "expand": True, "fatrop.mu_init":1e-1})
opti.to_function("opti_func", [], [opti.x]).generate('casadi_generated.c', {"with_header": True})
# opti.solve()