import numpy as np
import casadi as cs

def ldl_solve(A :cs.MX, b:cs.MX):
    def lin_solve_func(A:cs.MX, b:cs.MX):
        def subs_backward(U :cs.SX, y:cs.SX):
            m = U.shape[0]
            n = U.shape[1]
            x = cs.SX.zeros(n)
            x[:] = y[:]
            assert(m ==n)
            assert(y.shape[0] == m)
            for i in range(n-1, -1, -1):
                x[i] = x[i]/U[i,i]
                x[:i] -= U[:i,i]*x[i]
            return x

        def subs_forward(L :cs.SX, y:cs.SX):
            m = L.shape[0]
            n = L.shape[1]
            x = cs.SX.zeros(n)
            x[:] = y[:]
            assert(m ==n)
            assert(y.shape[0] == m)
            for i in range(n):
                x[i] = x[i]/L[i,i]
                x[i+1:] -= L[i+1:,i]*x[i]
            return x
        m = A.shape[0]
        n = A.shape[1]
        assert(m == n)
        assert(b.shape[0] == m)
        A_sym = cs.SX.sym('A', m, n)
        y_sym = cs.SX.sym('y', m)
        L = cs.chol(A_sym).T
        x = cs.SX.zeros(n)
        x[:] = y_sym[:]
        y =subs_backward(L.T, subs_forward(L, x))
        f = cs.Function('linear_solve_cust', [A_sym, y_sym], [y])
        return f
    return lin_solve_func(A, b)(A, b)

def get_ddq(q, dq, T, V, W):
    #  helper functions
    def time_der(expr):
        return [cs.jacobian(expr, q), cs.jacobian(expr, dq)] #  @ [dq, ddq]
    L = T - V + W
    dt_L_dq = time_der(cs.gradient(L, dq))
    A = dt_L_dq[1]
    b = cs.gradient(L, q) - dt_L_dq[0]@dq
    ddq = ldl_solve(A, b)
    return ddq