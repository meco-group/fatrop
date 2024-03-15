import numpy as np
import casadi as cs


def ldl_fact(A):
    def lin_fact_func(A):
        m = A.shape[0]
        A_sym = cs.SX.sym('A', m, m)
        L = cs.chol(A_sym).T
        f = cs.Function('ldl_fact_cust', [A_sym], [L])
        return f
    m = A.shape[0]
    assert(m == A.shape[1])
    return lin_fact_func(A)(A)


def ldl_solve(L:cs.MX, b:cs.MX):
    def subs_backward(U :cs.SX, y:cs.SX):
        m = U.shape[0]
        n = y.shape[1]
        assert(m ==U.shape[1])
        assert(y.shape[0] == m)
        x = cs.MX.zeros(m, n)
        x[:, :] = y[:, :]
        for i in range(m-1, -1, -1):
            x[i, :] = x[i, :]/U[i,i]
            x[:i, :] -= U[:i,i]@x[i, :]
        return x

    def subs_forward(L :cs.SX, y:cs.SX):
        m = L.shape[0]
        n = y.shape[1]
        assert(m == L.shape[1])
        assert(y.shape[0] == m)
        x = cs.MX.zeros(m, n)
        x[:, :] = y[:, :]
        for i in range(m):
            x[i, :] = x[i, :]/L[i,i]
            x[i+1:, :] -= L[i+1:,i]@x[i, :]
        return x
    m = L.shape[0]
    n = b.shape[1]
    assert(m == L.shape[1])
    assert(b.shape[0] == m)
    x = cs.MX.zeros(m, n)
    x[:,:] = b[:,:]
    y = subs_backward(L.T, subs_forward(L, x))
    return y
def symmeterize(A):
    m = A.shape[0]
    n = A.shape[1]
    assert(m == n)
    for i in range(m):
        for j in range(i-1):
            A[i,j] = A[j,i]
    return A


def get_ddq(q, dq, T, V, W):
    print("warning: this Lagrangian dynamics implementation is not very efficient (slow function evaluation), it is only for illustrative purposes. Use a more dynamics efficient method for real problems.")
    #  helper functions
    def time_der(expr):
        return [cs.jacobian(expr, q), cs.jacobian(expr, dq)] #  @ [dq, ddq]
    L = T - V + W
    dt_L_dq = time_der(cs.gradient(L, dq))
    A = dt_L_dq[1]
    b = cs.gradient(L, q) - dt_L_dq[0]@dq
    ddq = ldl_solve(ldl_fact(A), b)
    return ddq

def get_ddq_jac(q, dq, u, T, V, W):
    print("warning: this Lagrangian dynamics implementation is not very efficient (slow function evaluation), it is only for illustrative purposes. Use a more dynamics efficient method for real problems.")
    def time_der(expr):
        return [cs.jacobian(expr, q), cs.jacobian(expr, dq)] #  @ [dq, ddq]
    L = T - V + W
    dt_L_dq = time_der(cs.gradient(L, dq))
    A = symmeterize(dt_L_dq[1])
    b = cs.gradient(L, q) - dt_L_dq[0]@dq
    L_fact = ldl_fact(A)
    ddq = ldl_solve(L_fact, b)
    uqdq = cs.vertcat(u, q, dq)
    ddq_dum = cs.MX.sym('ddq', ddq.shape[0], ddq.shape[1])
    jac =cs.substitute(ldl_solve(L_fact, cs.jacobian(b - A@ddq_dum, uqdq)), ddq_dum, ddq)

    # g = A ddq - b = 0
    # x = uqdq    y = ddq
    # dydx = jac
    lam = cs.MX.sym('lam', ddq.shape[0])
    f = ddq_dum.T @ lam
    mu = ldl_solve(L_fact, cs.gradient(f, ddq_dum))
    z = cs.vertcat(uqdq, ddq_dum)
    mu_dummy = cs.MX.sym('mu', mu.shape[0], mu.shape[1])
    h_partial= cs.gradient(f, z) - (cs.jacobian(A@ddq_dum - b, z)).T@mu_dummy
    hess_full = cs.jacobian(h_partial, z)
    nx = uqdq.shape[0]
    hess_full, h = cs.substitute([hess_full, h_partial[:nx]], [mu_dummy, ddq_dum], [mu, ddq])
    hess = symmeterize(cs.horzcat(cs.MX.eye(nx), jac.T) @ hess_full @ cs.vertcat(cs.MX.eye(nx), jac))
    return ddq, jac, (hess, h, lam)