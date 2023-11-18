import casadi as cs
def transf(theta, p):
    return cs.vertcat(
        cs.horzcat(cs.cos(theta), -cs.sin(theta), p[0]),
                     cs.horzcat(cs.sin(theta), cs.cos(theta), p[1]),
                     cs.horzcat(0.0, 0.0, 1.))
