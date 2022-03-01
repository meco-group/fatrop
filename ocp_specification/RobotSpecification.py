
from OCPSpecification import *
import FSDynamics
from casadi import *
from urdf2casadi import urdfparser as u2c
ind_lower = np.array([1, 2, 5])
class RobotSpecification(OCPSpecificationInterface):
    def __init__(self, w_pos=1, w_rot=1, w_invars=(10**-3)*np.array([1.0, 1.0, 1.0])):
        root_link = "panda_link0"
        end_link = "panda_hand"
        robot_parser = u2c.URDFparser()
        robot_parser.from_file("./panda_arm_model.urdf")
        # Also supports .from_server for ros parameter server, or .from_string if you have the URDF as a string.
        fk_dict = robot_parser.get_forward_kinematics(root_link, end_link)
        self.n_joints = robot_parser.get_n_joints(root_link, end_link)
        # should give ['q', 'upper', 'lower', 'dual_quaternion_fk', 'joint_names', 'T_fk', 'joint_list', 'quaternion_fk']
        self.forward_kinematics = fk_dict["T_fk"]
        self.lower = np.array(fk_dict["lower"])
        self.upper = np.array(fk_dict["upper"])
        super().__init__()
        return

    def SetProblemDimensions(self):
        self.nx = self.n_joints 
        self.nu = self.n_joints 
        self.ngI = self.n_joints 
        self.ngF =  3 
        self.ngIneq = self.n_joints 
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
        return 0.5*(self.lower + self.upper) - xk

    def EqConstrFinal(self, xK, stage_params, global_params):
        return self.forward_kinematics(xK)[:3,3] - global_params

    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        return [self.lower, xk, self.upper]
