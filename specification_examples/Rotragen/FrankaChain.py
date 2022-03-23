import numpy as np
import casadi as cs
from urdf2casadi import urdfparser as u2c
from Capsule import Capsule, cs2np
import urdf_parser_py.urdf as urdf
from typing import List


class FrankaChain:
    # link_list should have exactly one joint in between and be in the right order
    def __init__(self, root_link, link_list, end_effector, urdf_file):
        self.root_link_ = root_link
        self.link_list_ = link_list
        self.robot_parser = u2c.URDFparser()
        self.robot_parser.from_file(urdf_file)
        self.n_joints = len(link_list)
        self.n_end_efs = len(end_effector)
        self.fk_dict = self.robot_parser.get_forward_kinematics(root_link, end_effector[0])
        self.joint_lower = np.array(self.fk_dict["lower"])
        self.joint_upper = np.array(self.fk_dict["upper"])
        self.end_effector_list_ = end_effector
        f = open(urdf_file, 'r')
        self.robot = urdf.Robot.from_xml_string(f.read())
        return

    def set_up_expressions(self, joint_states_sym):
        # set up link_expressions
        # order of this list is same as link_list
        # self.joint_angles_sym = cs.SX.sym('joinst_states', self.n_joints)
        self.frame_expr_list = [self.robot_parser.get_forward_kinematics(
            self.root_link_, self.link_list_[0])["T_fk"](joint_states_sym[0])]
        self.frame_expr_dict = {}
        for i in range(1, self.n_joints):
            self.frame_expr_list.append(self.frame_expr_list[-1]@self.robot_parser.get_forward_kinematics(
                self.link_list_[i-1], self.link_list_[i])["T_fk"](joint_states_sym[i]))
        for i in range(self.n_end_efs):
            self.frame_expr_list.append(self.frame_expr_list[self.n_joints-1]@self.robot_parser.get_forward_kinematics(
            self.link_list_[self.n_joints-1], self.end_effector_list_[i])["T_fk"]()["o0"])
        for i in range(0, self.n_joints):
            self.frame_expr_dict[self.link_list_[i]] = self.frame_expr_list[i]
        for i in range(0, self.n_end_efs):
            self.frame_expr_dict[self.end_effector_list_[i]] = self.frame_expr_list[self.n_joints+i]
        self.capsulelist: List[Capsule] = []
        combined = self.link_list_ + self.end_effector_list_
        for link in self.robot.links:
            if link.name in combined:
                link_frame = self.GetFrameExpression(link.name)
                for coll in link.collisions:
                    if isinstance(coll.geometry, urdf.Cylinder):
                        cyl: urdf.Cylinder = coll.geometry
                        origin: urdf.Pose = coll.origin
                        print('adding collision for link ', link.name, ' xyz ', origin.xyz, ' rpy ', origin.rpy, ' r ', cyl.radius, ' L ', cyl.length)
                        self.capsulelist.append(
                            Capsule(link_frame, origin.xyz, origin.rpy, cyl.radius, cyl.length))

    def GetFrameIndex(self, frame_name):
        return self.link_list_.index(frame_name)

    def GetFrameExpression(self, frame_name):
        return self.frame_expr_dict[frame_name]

    def GetAllCapsuleConstraints(self, point, distance):
        res = cs.MX.zeros(len(self.capsulelist)) 
        for capsi in range(len(self.capsulelist)):
            caps = self.capsulelist[capsi]
            res[capsi] = caps.point_distance_constraint(point, distance)
        return res

if __name__ == "__main__":
    frankachain = FrankaChain('world', ['panda_link1', 'panda_link2', 'panda_link3', 'panda_link4',
                              'panda_link5', 'panda_link6', 'panda_link7'], ['TCP_frame', 'panda_link8', 'panda_hand'], 'panda_arm_model.urdf')

    obstablepos = np.array([0.6, 0.00,1.0]) 
    obstaclerad = 0.20
    frankachain.set_up_expressions(np.array([.0, 0.0,0.0,-1.0,0.0,2.5,00.0]))
    print(cs2np(frankachain.GetAllCapsuleConstraints(obstablepos, obstaclerad)))
    print(cs2np(frankachain.GetFrameExpression('TCP_frame')))

    # ####### solve an inverse kinematics problem
    # opti = cs.Opti()
    # joint_states = opti.variable(7)
    # obstacle = np.array([0.5, 0.0, 0.8])
    # radius = 0.1
    # target = np.array([0.8, 0.0, 0.5])
    # frankachain.set_up_expressions(joint_states)
    # opti.subject_to(frankachain.joint_lower<(joint_states<frankachain.joint_upper))
    # opti.subject_to(frankachain.GetAllCapsuleConstraints(obstacle, radius)>0)
    # alpha :float = 0.5
    # opti.set_initial(joint_states, alpha*frankachain.joint_upper + (1-alpha)*frankachain.joint_lower)
    # opti.subject_to(frankachain.GetFrameExpression('TCP_frame')[:3,3]-target == 0)
    # opti.subject_to(frankachain.GetFrameExpression('TCP_frame')[1:3,0] == 0)
    # opti.minimize(0.5*joint_states.T@joint_states)
    # opti.solver('ipopt')
    # res = opti.solve()
    # print(res.value(frankachain.GetAllCapsuleConstraints(obstacle, radius)))