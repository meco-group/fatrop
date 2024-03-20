import fatropy.spectool as sp
import casadi as cs
import numpy as np
dt = 0.01

ocp = sp.Ocp()
x = ocp.state()
v = ocp.state()
F = ocp.control()
x0 = ocp.parameter()
v0 = ocp.parameter()
stage = ocp.new_stage(50)
stage.set_next(v, v + dt * F)
stage.set_next(x, x + dt * v)
stage.add_objective(cs.sumsqr(F), sp.t0, sp.mid, sp.tf)
ocp.at_t0().subject_to(x == x0)
ocp.at_t0().subject_to(v == v0) 
ocp.solver("fatrop")
func = ocp.to_function("double_integrator", [x0, v0], [ocp.at_t0(F)])
func(0.1, 0.1)
