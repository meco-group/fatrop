import casadi as cs
import fatropy.spectool as sp

def gauss_newton_hessian(expr:cs.MX,vars:cs.MX):
    """
    Compute the Gauss-Newton approximation of the Hessian of the expression `expr` with respect to the decision variables.
    This is the Hessian that omits the second order terms of the Taylor expansion of the expression.
    """
    lam = cs.MX.sym("lam", expr.shape[0])
    return sp.Hessian(cs.MX.zeros(vars.shape[0], vars.shape[0]), cs.gradient(expr.T@lam, vars), vars, lam)