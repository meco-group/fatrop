import casadi as cs
import numpy as np
import lagrangian_dynamics
import robot2d as r2d


m_leg_lower = 1.
l_leg_lower = .25
m_leg_upper = 1.
l_lleg_upper = .25
m_body = 5.
l_body = .4
m_foot = .5
l_foot = 0.05

def ode_left_contact(q, dq, tau, x_foot_left):
    theta_foot_left, dtheta_foot_left = q[0], dq[0]
    theta_foot_right, dtheta_foot_right = q[1], dq[1]
    theta_knee_left, dtheta_knee_left = q[2], dq[2]
    theta_knee_right, dtheta_knee_right = q[3], dq[3]
    theta_hip_left, dtheta_hip_left = q[4], dq[4]
    theta_hip_right, dtheta_hip_right = q[5], dq[5]
    # built up the mechanism
    foot_contact_p = r2d.transform2d(x_foot_left, 0, 0)
    foot_left_joint = r2d.revolute(foot_contact_p.transform(r2d.transform2d(0, 0, np.pi/2)), theta_foot_left, dtheta_foot_left)
    lower_leg_left = r2d.link(foot_left_joint, m_leg_upper, l_lleg_upper)
    knee_left_joint = r2d.revolute(lower_leg_left.right, theta_knee_left, dtheta_knee_left)
    upper_leg_left = r2d.link(knee_left_joint, m_leg_upper, l_lleg_upper)
    hip_left_joint = r2d.revolute(upper_leg_left.right, theta_hip_left, dtheta_hip_left)
    body = r2d.link(hip_left_joint, m_body, l_body)
    hip_right_joint = r2d.revolute(body.left.transform(r2d.transform2d(0, 0, np.pi)), theta_hip_right, dtheta_hip_right)
    upper_leg_right = r2d.link(hip_right_joint, m_leg_upper, l_lleg_upper)
    knee_right_joint = r2d.revolute(upper_leg_right.right, theta_knee_right, dtheta_knee_right)
    lower_leg_right = r2d.link(knee_right_joint, m_leg_upper, l_lleg_upper)
    foot_right_joint = r2d.revolute(lower_leg_right.right.transform(r2d.transform2d(0, 0, np.pi/2)), theta_foot_right, dtheta_foot_right)
    foot_right = r2d.link(foot_right_joint, m_foot, l_foot)
    # collect all the links
    links = [lower_leg_left, upper_leg_left, body, upper_leg_right, lower_leg_right, foot_right]
    # W = F_left * lower_leg_left.left.x + M_left*lower_leg_left.left.theta
    W = tau.T@q
    # built up the mechanism
    mechanism = r2d.mechanism(links)
    ddq, _, _ = mechanism.get_dynamics(W, q, dq, tau)
    return ddq, mechanism, links

def ode_right_contact(q, dq, tau, x_foot_right):
    theta_foot_left, dtheta_foot_left = q[0], dq[0]
    theta_foot_right, dtheta_foot_right = q[1], dq[1]
    theta_knee_left, dtheta_knee_left = q[2], dq[2]
    theta_knee_right, dtheta_knee_right = q[3], dq[3]
    theta_hip_left, dtheta_hip_left = q[4], dq[4]
    theta_hip_right, dtheta_hip_right = q[5], dq[5]
    # built up the mechanism
    foot_contact_p = r2d.transform2d(x_foot_right, 0, 0.)
    foot_right_joint = r2d.revolute(foot_contact_p.transform(r2d.transform2d(0, 0, np.pi/2)), -theta_foot_right, -dtheta_foot_right)
    lower_leg_right = r2d.link(foot_right_joint, m_leg_upper, l_lleg_upper)
    knee_right_joint = r2d.revolute(lower_leg_right.right, -theta_knee_right, -dtheta_knee_right)
    upper_leg_right = r2d.link(knee_right_joint, m_leg_upper, l_lleg_upper)
    hip_right_joint = r2d.revolute(upper_leg_right.right, -theta_hip_right, -dtheta_hip_right)
    body = r2d.link(hip_right_joint, m_body, l_body)
    hip_left_joint = r2d.revolute(body.left.transform(r2d.transform2d(0, 0, -np.pi)), -theta_hip_left, -dtheta_hip_left)
    upper_leg_left = r2d.link(hip_left_joint, m_leg_upper, l_lleg_upper)
    knee_left_joint = r2d.revolute(upper_leg_left.right, -theta_knee_left, -dtheta_knee_left)
    lower_leg_left = r2d.link(knee_left_joint, m_leg_upper, l_lleg_upper)
    foot_left_joint = r2d.revolute(lower_leg_left.right.transform(r2d.transform2d(0, 0, np.pi/2)), -theta_foot_left, -dtheta_foot_left)
    foot_left = r2d.link(foot_left_joint, m_foot, l_foot)
    links = [lower_leg_left, upper_leg_left, body, upper_leg_right, lower_leg_right, foot_left]
    # W = F_left * lower_leg_left.left.x + M_left*lower_leg_left.left.theta
    W = tau.T@q
    # built up the mechanism
    mechanism = r2d.mechanism(links)
    ddq, _, _ = mechanism.get_dynamics(W, q, dq, tau)
    return ddq, mechanism, links


def ode_left_contact_arms(q, dq, tau, x_foot_left):
    theta_foot_left, dtheta_foot_left = q[0], dq[0]
    theta_foot_right, dtheta_foot_right = q[1], dq[1]
    theta_knee_left, dtheta_knee_left = q[2], dq[2]
    theta_knee_right, dtheta_knee_right = q[3], dq[3]
    theta_hip_left, dtheta_hip_left = q[4], dq[4]
    theta_hip_right, dtheta_hip_right = q[5], dq[5]
    shoulder_left, dshoulder_left = q[6], dq[6]
    shoulder_right, dshoulder_right = q[7], dq[7]
    elbow_left, delbow_left = q[8], dq[8]
    elbow_right, delbow_right = q[9], dq[9]
    # built up the mechanism
    foot_contact_p = r2d.transform2d(x_foot_left, 0, 0)
    foot_left_joint = r2d.revolute(foot_contact_p.transform(r2d.transform2d(0, 0, np.pi/2)), theta_foot_left, dtheta_foot_left)
    lower_leg_left = r2d.link(foot_left_joint, m_leg_upper, l_lleg_upper)
    knee_left_joint = r2d.revolute(lower_leg_left.right, theta_knee_left, dtheta_knee_left)
    upper_leg_left = r2d.link(knee_left_joint, m_leg_upper, l_lleg_upper)
    hip_left_joint = r2d.revolute(upper_leg_left.right, theta_hip_left, dtheta_hip_left)
    body = r2d.link(hip_left_joint, m_body, l_body)
    hip_right_joint = r2d.revolute(body.left.transform(r2d.transform2d(0, 0, np.pi)), theta_hip_right, dtheta_hip_right)
    upper_leg_right = r2d.link(hip_right_joint, m_leg_upper, l_lleg_upper)
    knee_right_joint = r2d.revolute(upper_leg_right.right, theta_knee_right, dtheta_knee_right)
    lower_leg_right = r2d.link(knee_right_joint, m_leg_upper, l_lleg_upper)
    foot_right_joint = r2d.revolute(lower_leg_right.right.transform(r2d.transform2d(0, 0, np.pi/2)), theta_foot_right, dtheta_foot_right)
    foot_right = r2d.link(foot_right_joint, m_foot, l_foot)
    arm_left_joint = r2d.revolute(body.right.transform(r2d.transform2d(0, 0, np.pi)), shoulder_left, dshoulder_left)
    arm_left = r2d.link(arm_left_joint, m_leg_upper, l_lleg_upper)
    arm_right_joint = r2d.revolute(body.right.transform(r2d.transform2d(0, 0, np.pi)), shoulder_right, dshoulder_right)
    arm_right = r2d.link(arm_right_joint, m_leg_upper, l_lleg_upper)
    arm2_left_joint = r2d.revolute(arm_left.right, elbow_left, delbow_left)
    arm2_left = r2d.link(arm2_left_joint, m_leg_upper, l_lleg_upper)
    arm2_right_joint = r2d.revolute(arm_right.right, elbow_right, delbow_right)
    arm2_right = r2d.link(arm2_right_joint, m_leg_upper, l_lleg_upper)
    # collect all the links
    links = [arm_left, arm_right, arm2_left, arm2_right, lower_leg_left, upper_leg_left, body, upper_leg_right, lower_leg_right, foot_right]
    # W = F_left * lower_leg_left.left.x + M_left*lower_leg_left.left.theta
    W = tau.T@q
    # built up the mechanism
    mechanism = r2d.mechanism(links)
    ddq, _, _ = mechanism.get_dynamics(W, q, dq, tau)
    return ddq, mechanism, links

def ode_right_contact_arms(q, dq, tau, x_foot_right):
    theta_foot_left, dtheta_foot_left = q[0], dq[0]
    theta_foot_right, dtheta_foot_right = q[1], dq[1]
    theta_knee_left, dtheta_knee_left = q[2], dq[2]
    theta_knee_right, dtheta_knee_right = q[3], dq[3]
    theta_hip_left, dtheta_hip_left = q[4], dq[4]
    theta_hip_right, dtheta_hip_right = q[5], dq[5]
    shoulder_left, dshoulder_left = q[6], dq[6]
    shoulder_right, dshoulder_right = q[7], dq[7]
    elbow_left, delbow_left = q[8], dq[8]
    elbow_right, delbow_right = q[9], dq[9]
    # built up the mechanism
    foot_contact_p = r2d.transform2d(x_foot_right, 0, 0.)
    foot_right_joint = r2d.revolute(foot_contact_p.transform(r2d.transform2d(0, 0, np.pi/2)), -theta_foot_right, -dtheta_foot_right)
    lower_leg_right = r2d.link(foot_right_joint, m_leg_upper, l_lleg_upper)
    knee_right_joint = r2d.revolute(lower_leg_right.right, -theta_knee_right, -dtheta_knee_right)
    upper_leg_right = r2d.link(knee_right_joint, m_leg_upper, l_lleg_upper)
    hip_right_joint = r2d.revolute(upper_leg_right.right, -theta_hip_right, -dtheta_hip_right)
    body = r2d.link(hip_right_joint, m_body, l_body)
    hip_left_joint = r2d.revolute(body.left.transform(r2d.transform2d(0, 0, -np.pi)), -theta_hip_left, -dtheta_hip_left)
    upper_leg_left = r2d.link(hip_left_joint, m_leg_upper, l_lleg_upper)
    knee_left_joint = r2d.revolute(upper_leg_left.right, -theta_knee_left, -dtheta_knee_left)
    lower_leg_left = r2d.link(knee_left_joint, m_leg_upper, l_lleg_upper)
    foot_left_joint = r2d.revolute(lower_leg_left.right.transform(r2d.transform2d(0, 0, np.pi/2)), -theta_foot_left, -dtheta_foot_left)
    foot_left = r2d.link(foot_left_joint, m_foot, l_foot)
    arm_left_joint = r2d.revolute(body.right.transform(r2d.transform2d(0, 0, np.pi)), shoulder_left, dshoulder_left)
    arm_left = r2d.link(arm_left_joint, m_leg_upper, l_lleg_upper)
    arm_right_joint = r2d.revolute(body.right.transform(r2d.transform2d(0, 0, np.pi)), shoulder_right, dshoulder_right)
    arm_right = r2d.link(arm_right_joint, m_leg_upper, l_lleg_upper)
    arm2_left_joint = r2d.revolute(arm_left.right, elbow_left, delbow_left)
    arm2_left = r2d.link(arm2_left_joint, m_leg_upper, l_lleg_upper)
    arm2_right_joint = r2d.revolute(arm_right.right, elbow_right, delbow_right)
    arm2_right = r2d.link(arm2_right_joint, m_leg_upper, l_lleg_upper)
    links = [arm_left, arm_right, arm2_left, arm2_right, lower_leg_left, upper_leg_left, body, upper_leg_right, lower_leg_right, foot_left]
    # W = F_left * lower_leg_left.left.x + M_left*lower_leg_left.left.theta
    W = tau.T@q
    # built up the mechanism
    mechanism = r2d.mechanism(links)
    ddq, _, _ = mechanism.get_dynamics(W, q, dq, tau)
    return ddq, mechanism, links