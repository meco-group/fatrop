import walker2d_dynamics
import casadi as cs
import numpy as np
import fatropy.spectool as sp
import robot2d as r2d

N = 10
N_steps = 4
T = .5
ocp = sp.Ocp()

q = ocp.state(6)
qd = ocp.state(6)
tau = ocp.control(6)

q0 = cs.DM([0, 0, 0, 0, 0, 0.]) # foot_left, foot_right, knee_left, knee_right, hip_left, hip_right
qd0 = cs.DM([0, 0, 0, 0, 0, 0.])
qf = cs.DM([0, 0., 0, 0, 0, 0.])
qdf = cs.DM([0, 0, 0, 0, 0, 0.])

all_steps = []
all_mechanisms = []

for step in range(N_steps):
    ddq_left, mechanism_left, links_left = walker2d_dynamics.ode_left_contact(q, qd, tau, 2*step*0.25-1.)
    ddq_right, mechanism_right, links_right = walker2d_dynamics.ode_right_contact(q, qd, tau, (2*step+1)*0.25-1.)
    all_mechanisms += (mechanism_left, mechanism_right)
    stage_left = ocp.new_stage(N)
    stage_right = ocp.new_stage(N)
    all_steps.append((stage_left, stage_right))
    # Define the dynamics

    stage_left.set_next(q, q + qd * T / N)
    stage_left.set_next(qd, qd + ddq_left * T / N)

    stage_right.set_next(q, q + qd * T / N)
    stage_right.set_next(qd, qd + ddq_right * T / N)

    # transition dynamics
    stage_left.at_tf().set_next(q, q)
    stage_left.at_tf().set_next(qd, qd)
    if step != N_steps - 1:
        stage_right.at_tf().set_next(q, q)
        stage_right.at_tf().set_next(qd, qd)

    # transition constraints
    foot_right = links_left[-1]
    stage_left.at_tf().subject_to(foot_right.left.origin() == [(2*step+1)*0.25-1., 0.0])
    stage_left.at_tf().subject_to(cs.sin(foot_right.left.theta) == 0.) # theta = 0. + k * pi
    if step != N_steps - 1:
        foot_left = links_right[-1]
        stage_right.at_tf().subject_to(foot_left.left.origin() == [(2*step+2)*0.25-1., 0.0])
        stage_right.at_tf().subject_to(cs.sin(foot_left.left.theta) == 0.) # theta = 0. + k * pi

    # Define the cost function
    for stage in [stage_left, stage_right]:
        stage.add_objective(cs.sumsqr(qd), sp.t0, sp.mid, sp.tf)
        stage.add_objective(10*cs.sumsqr(tau), sp.t0, sp.mid)
    # no over-bending of knee
    for stage in [stage_left, stage_right]:
        stage.subject_to(0. <= q[2], sp.mid)
        stage.subject_to(0. >= q[3], sp.mid)
    # no torque in wrong direction from feet
    stage_left.subject_to(tau[0]<0., sp.mid)
    stage_right.subject_to(tau[1]>0., sp.mid)


# Define the constraints
ocp.at_t0().subject_to(q == q0)
ocp.at_t0().subject_to(qd == qd0)
ocp.at_tf().subject_to(q == qf)
ocp.at_tf().subject_to(qd == qdf)


ocp.solver('fatrop', {"post_expand":True, "jit":False}, {})

fun = ocp.to_function("fun", [], [ocp.sample(q)[1], stage_left.at_tf(foot_right.left.origin())])

sol = fun()
q_sol = sol['o0']
r2d.animate([(all_mechanisms[i], q, q_sol[:, i*(N+1):(i+1)*(N+1)-1]) for i in range(2*N_steps)])