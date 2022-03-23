from tkinter import Y
import numpy as np
import numpy.linalg
import casadi as cs
import urdf2casadi
def cs2np(asd):
    return cs.Function("temp",[],[asd])()["o0"].toarray()
def inv_transfmatrix(transfm):
    res = cs.MX.zeros(4,4)
    R = transfm[:3,:3]
    p = transfm[:3, 3]
    res[:3,:3] = R.T
    res[:3,3] = - R.T@p
    res[3,3] = 1.0
    return res
class Capsule:
    def __init__(self, frame: cs.MX, origin: np.array, rpy: np.array, radius: float, length: float):
        self.frame = frame
        self.rpy = rpy
        self.origin_ = origin
        self.radius = radius 
        self.length = length
        roll: float = rpy[0]
        pitch: float = rpy[1]
        yaw: float = rpy[2]
        x: float = origin[0]
        y: float = origin[1]
        z: float = origin[2]
        self.T_world_capsule = frame @ cs2np(urdf2casadi.geometry.transformation_matrix.full_symbolic([x,y,z], [roll, pitch, yaw]))  ## this is a casadi expression
        self.T_capsule_world = inv_transfmatrix(self.T_world_capsule)
    def point_distance_constraint_old(self, point, min_dist):
        ## the value of this expression should be larger than zero
        ## point coordinates in local frame
        point_local = (self.T_capsule_world @ cs.vertcat(point, np.array([1.0])))[:3]
        x_l = point_local[0]
        y_l = point_local[1]
        z_l = cs.if_else(point_local[2]<0, -1, 1)*point_local[2]
        ## distance to cylinder vs distance to ball
        return cs.if_else(z_l<self.length/2, x_l*x_l+ y_l*y_l - (min_dist+ self.radius)**2, x_l**2+y_l**2+(z_l-self.length/2)**2 - (min_dist+ self.radius)**2)
    def point_distance_constraint(self, point, min_dist):
        ## the value of this expression should be larger than zero
        ## point coordinates in local frame
        point_local = (self.T_capsule_world @ cs.vertcat(point, np.array([1.0])))[:3]
        x_l = point_local[0]
        y_l = point_local[1]
        z_l = cs.if_else(point_local[2]<0, -1, 1)*point_local[2]
        ## distance to cylinder vs distance to ball
        return cs.if_else(z_l<self.length/2, cs.sqrt(x_l*x_l+ y_l*y_l) - min_dist- self.radius, cs.sqrt(x_l*x_l+y_l*y_l+(z_l-self.length/2)**2)- (min_dist+ self.radius))
    def point_distance_constraint_ball(self, point, min_dist):
        ## the value of this expression should be larger than zero
        ## point coordinates in local frame
        origin_caps = self.T_world_capsule[:3,3]
        delta = point - origin_caps
        return (cs.sum1(delta**2)- (min_dist+self.length+self.radius)**2)
    def point_distance_constraint_old(self, point, min_dist):
        ## the value of this expression should be larger than zero
        ## point coordinates in local frame
        origin_caps = self.T_world_capsule[:3,3]
        delta = point - origin_caps
        projected_length = delta.T @ self.T_world_capsule[:3,2]
        projected_length_abs = cs.if_else(projected_length<0, -projected_length, projected_length)
        delta_or = delta - (projected_length * self.T_world_capsule[:3,2])
        delta_origin_ball = self.length/2 *self.T_world_capsule[:3,2]
        or1 = origin_caps+delta_origin_ball
        or2 = origin_caps-delta_origin_ball
        mindistr = (min_dist+self.radius)**2
        length1 = cs.sum1((point-or1)**2)- mindistr
        length2 = cs.sum1((point-or2)**2)- mindistr
        return cs.if_else(projected_length_abs< self.length/2,cs.sum1(delta_or**2) - mindistr, cs.if_else(length1<length2, length1, length2))



