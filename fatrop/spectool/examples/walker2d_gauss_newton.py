import walker2d_dynamics
import casadi as cs
import numpy as np
import fatropy.spectool as sp
import robot2d as r2d
from gauss_newton_hessian import gauss_newton_hessian

N = 50
T = .5
ocp = sp.Ocp()

stage_left = ocp.new_stage(N)
stage_right = ocp.new_stage(N)
q = ocp.state(6)
qd = ocp.state(6)
tau = ocp.control(6)
ux = cs.vertcat(tau, q, qd)
x = cs.vertcat(q, qd)

q0 = cs.DM([0, 0, 0, 0, 0, 0.]) # foot_left, foot_right, knee_left, knee_right, hip_left, hip_right
qd0 = cs.DM([0, 0, 0, 0, 0, 0.])
qf = cs.DM([0, 0., 0, 0, 0, 0.])
qdf = cs.DM([0, 0, 0, 0, 0, 0.])


# Define the dynamics
ddq_left, mechanism_left, links_left = walker2d_dynamics.ode_left_contact(q, qd, tau, 0.0)
ddq_right, mechanism_right, links_right = walker2d_dynamics.ode_right_contact(q, qd, tau, 0.25)

stage_left.set_next(q, q + qd * T / N)
qd_next = qd + ddq_left * T / N
stage_left.set_next(qd, qd_next, hessian = gauss_newton_hessian(qd_next, ux))

stage_right.set_next(q, q + qd * T / N)
qd_next = qd + ddq_right * T / N
stage_right.set_next(qd, qd_next, hessian = gauss_newton_hessian(qd_next, ux))


# Define the constraints
ocp.at_t0().subject_to(q == q0)
ocp.at_t0().subject_to(qd == qd0)
ocp.at_tf().subject_to(q == qf)
ocp.at_tf().subject_to(qd == qdf)
# no over-bending of knees
for stage in [stage_left, stage_right]:
    stage.subject_to(0. <= q[2], sp.mid)
    stage.subject_to(0. >= q[3], sp.mid)

# transition dynamics
stage_left.at_tf().set_next(q, q)
stage_left.at_tf().set_next(qd, qd)

# transition constraints
foot_right = links_left[-1]
expr = foot_right.left.origin() -[0.25, 0.0]
stage_left.at_tf().subject_to(expr == 0., hessian = gauss_newton_hessian(expr, x))
expr = cs.sin(foot_right.left.theta) # theta = 0. + k * pi
stage_left.at_tf().subject_to(expr == 0., hessian = gauss_newton_hessian(expr, x)) 

# Define the cost function
for stage in [stage_left, stage_right]:
    stage.add_objective(cs.sumsqr(qd), sp.t0, sp.mid, sp.tf)
    stage.add_objective(1e-1*cs.sumsqr(tau), sp.t0, sp.mid)

ocp.solver('fatrop', {"post_expand":True, "jit":False}, {})

fun = ocp.to_function("fun", [], [ocp.sample(q)[1], stage_left.at_tf(foot_right.left.origin())])

sol = fun()
q_sol = sol['o0']
r2d.animate([(mechanism_left, q, q_sol[:, :N]), (mechanism_right, q, q_sol[:, N+2:])])