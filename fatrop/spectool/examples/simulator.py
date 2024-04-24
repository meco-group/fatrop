import casadi as cs
import numpy as np
from typing  import List

class simulator:
    def __init__(self, x: cs.MX, u :cs.MX, xp :cs.MX):
        self.x = x
        self.u = u
        self.xp = xp
        self.f = cs.Function('f', [x, u], [xp])

    def step(self, x, u):
        return self.f(x, u)

    def simulate(self, x0 : cs.DM , u: cs.DM, out: List[cs.MX]):
        N = u.shape[1]
        x_res = [x0]
        u = cs.horzsplit(u)
        for i in range(N):
            x_res.append(self.step(x_res[-1], u[i]))
        out_map = cs.Function('out_map', [self.x, self.u], out)
        u += [u[-1]] # to make the length of u equal to N
        return [cs.horzcat(*outi)  for outi in zip(*[out_map(x_res[i], u[i]) for i in range(N+1)])]