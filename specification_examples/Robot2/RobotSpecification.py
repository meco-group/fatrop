from numpy import Inf
import numpy.random
from FatropOCPSpecification import *
from casadi import *
from urdf2casadi import urdfparser as u2c

ind_lower = np.array([1, 2, 5])
n_points = 1
radius = np.array(n_points*[0.15, 0.12,0.1, 0.14,0.08, 0.14])
class RobotSpecification(OCPSpecificationInterface):
    def __init__(self, w_pos=1, w_rot=1, w_invars=(10**-3)*np.array([1.0, 1.0, 1.0])):
        root_link = "panda_link0"
        end_link = "panda_hand"
        self.robot_parser = u2c.URDFparser()
        self.robot_parser.from_file("./panda_arm_model.urdf")
        # Also supports .from_server for ros parameter server, or .from_string if you have the URDF as a string.
        self.fk_dict = self.robot_parser.get_forward_kinematics(root_link, end_link)
        self.n_joints = self.robot_parser.get_n_joints(root_link, end_link)
        print('n_joints')
        print(self.n_joints)
        # should give ['q', 'upper', 'lower', 'dual_quaternion_fk', 'joint_names', 'T_fk', 'joint_list', 'quaternion_fk']
        self.forward_kinematics = self.fk_dict["T_fk"]
        self.joint_lower = np.array(self.fk_dict["lower"])
        self.joint_upper = np.array(self.fk_dict["upper"])
        print("lower " , self.joint_lower)
        print("upper " , self.joint_upper)
        max_vel = pi # rad/sec
        self.lower = np.hstack((self.joint_lower, -max_vel*np.ones(self.n_joints), radius**2))
        self.upper = np.hstack((self.joint_upper, max_vel *np.ones(self.n_joints),np.array(n_points*6*[1e5])))
        # self.lowerF = np.hstack((self.joint_lower, radius**2))
        # self.upperF = np.hstack((self.joint_upper, np.array(n_points*6*[1e5])))
        self.lowerF = self.joint_lower 
        self.upperF = self.joint_upper
        super().__init__()
        return

    def SetProblemDimensions(self):
        self.nx = self.n_joints 
        self.nu = self.n_joints 
        self.ngI = self.n_joints 
        self.ngF =  3 
        self.ngIneq = self.n_joints+ self.n_joints + n_points*6
        # self.ngIneqF = self.n_joints + n_points*6
        self.ngIneqF = self.n_joints 
        self.n_stage_params = 1  # dt
        self.n_global_params = 3 # endpos 

    def Dynamics(self, uk, xk, stage_params, global_params):
        return xk + stage_params[0]*uk

    def StageCost(self, uk, xk, stage_params, global_params):
        ukxk = vertcat(uk, xk)
        H = np.eye(self.nx+self.nu)
        return 0.5 * (ukxk.T@H@ukxk)

    def StageCostFinal(self, xK, stage_params, global_params):
        H = np.eye(self.nx)
        return 0.5 * (xK.T@H@xK)
    
    def EqConstrInitial(self, uk, xk, stage_params, global_params):
        # return 0.5*(self.joint_lower + self.joint_upper) - xk
        return np.array(3*[0.0] + [-1.0] + 3*[0.0]) - xk

    def EqConstrFinal(self, xK, stage_params, global_params):
        return self.forward_kinematics(xK)[:3,3] - global_params

    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        root_link = "panda_link0"
        joint_expr0 = np.eye(4)
        joint_expr1 = joint_expr0 @ self.robot_parser.get_forward_kinematics("panda_link0", "panda_link1")["T_fk"](xk[0])
        joint_expr2 = joint_expr1 @ self.robot_parser.get_forward_kinematics("panda_link1", "panda_link2")["T_fk"](xk[1])
        joint_expr3 = joint_expr2 @ self.robot_parser.get_forward_kinematics("panda_link2", "panda_link3")["T_fk"](xk[2])
        joint_expr4 = joint_expr3 @ self.robot_parser.get_forward_kinematics("panda_link3", "panda_link4")["T_fk"](xk[3])
        joint_expr5 = joint_expr4 @ self.robot_parser.get_forward_kinematics("panda_link4", "panda_link5")["T_fk"](xk[4])
        joint_expr6 = joint_expr5 @ self.robot_parser.get_forward_kinematics("panda_link5", "panda_link6")["T_fk"](xk[5])
        joint_expr7 = joint_expr6 @ self.robot_parser.get_forward_kinematics("panda_link6", "panda_link7")["T_fk"](xk[6])
        # joint_expr8 = joint_expr7 @ self.robot_parser.get_forward_kinematics("panda_link7", "panda_hand")["T_fk"](0)
        joint_pos1 = joint_expr1[:3,3]
        joint_pos2 = joint_expr2[:3,3]
        joint_pos3 = joint_expr3[:3,3]
        joint_pos4 = joint_expr4[:3,3]
        joint_pos5 = joint_expr5[:3,3]
        joint_pos6 = joint_expr6[:3,3]
        joint_pos7 = joint_expr7[:3,3]
        # joint_pos8 = joint_expr8[:3,3]
        jointposlist = [joint_pos2, joint_pos3, joint_pos4, joint_pos5, joint_pos6, joint_pos7]
        distancelist = SX.zeros(6)
        distance_all = SX.zeros(6*n_points)
        for point in range(n_points):
            obstaclepos = np.array([0.5, 0.0,-10]) +0.00* np.random.rand(3)
            for i in range(6):
                distance_all[point*6+i] = sum1((jointposlist[i]-obstaclepos)**2)
            # distance_all[point*6:(point+1)*6] = distancelist[:]
        joint_vel = uk/global_params[0]
        
        return [self.lower, vertcat(xk, joint_vel, distance_all), self.upper]
    def FinalInequality(self, xk, stage_params, global_params):
        root_link = "panda_link0"
        joint_expr0 = np.eye(4)
        joint_expr1 = joint_expr0 @ self.robot_parser.get_forward_kinematics("panda_link0", "panda_link1")["T_fk"](xk[0])
        joint_expr2 = joint_expr1 @ self.robot_parser.get_forward_kinematics("panda_link1", "panda_link2")["T_fk"](xk[1])
        joint_expr3 = joint_expr2 @ self.robot_parser.get_forward_kinematics("panda_link2", "panda_link3")["T_fk"](xk[2])
        joint_expr4 = joint_expr3 @ self.robot_parser.get_forward_kinematics("panda_link3", "panda_link4")["T_fk"](xk[3])
        joint_expr5 = joint_expr4 @ self.robot_parser.get_forward_kinematics("panda_link4", "panda_link5")["T_fk"](xk[4])
        joint_expr6 = joint_expr5 @ self.robot_parser.get_forward_kinematics("panda_link5", "panda_link6")["T_fk"](xk[5])
        joint_expr7 = joint_expr6 @ self.robot_parser.get_forward_kinematics("panda_link6", "panda_link7")["T_fk"](xk[6])
        # joint_expr8 = joint_expr7 @ self.robot_parser.get_forward_kinematics("panda_link7", "panda_hand")["T_fk"](0)
        joint_pos1 = joint_expr1[:3,3]
        joint_pos2 = joint_expr2[:3,3]
        joint_pos3 = joint_expr3[:3,3]
        joint_pos4 = joint_expr4[:3,3]
        joint_pos5 = joint_expr5[:3,3]
        joint_pos6 = joint_expr6[:3,3]
        joint_pos7 = joint_expr7[:3,3]
        # joint_pos8 = joint_expr8[:3,3]
        jointposlist = [joint_pos2, joint_pos3, joint_pos4, joint_pos5, joint_pos6, joint_pos7]
        distance_all = SX.zeros(6)
        obstaclepos = np.array([0.5, 0.0,-10.0]) +0.00* np.random.rand(3)
        for i in range(6):
            distance_all[i] = sum1((jointposlist[i]-obstaclepos)**2)
        
        # return [self.lowerF, vertcat(xk, distance_all), self.upperF]
        return [self.lowerF[:7], vertcat(xk), self.upperF[:7]]