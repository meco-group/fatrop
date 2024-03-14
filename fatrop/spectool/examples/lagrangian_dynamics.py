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
    def time_der(expr):
        return [cs.jacobian(expr, q), cs.jacobian(expr, dq)] #  @ [dq, ddq]
    L = T - V + W
    dt_L_dq = time_der(cs.gradient(L, dq))
    A = dt_L_dq[1]
    b = cs.gradient(L, q) - dt_L_dq[0]@dq
    L_fact = ldl_fact(A)
    ddq = ldl_solve(L_fact, b)
    uqdq = cs.vertcat(u, q, dq)
    # ddq_dum = cs.MX.sym('ddq', ddq.shape[0], ddq.shape[1])
    # jac =cs.substitute(ldl_solve(L_fact, cs.jacobian(b - A@ddq_dum, uqdq)), ddq_dum, ddq)
    jac = cs.jacobian(ddq, uqdq)
    return ddq, jac