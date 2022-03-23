from cmath import inf
from turtle import distance
from numpy import Inf
import numpy.matlib
import numpy.random
from FatropOCPSpecification import *
from casadi import *
from urdf2casadi import urdfparser as u2c
from FrankaChain import FrankaChain

class RobSpecCaps(OCPSpecificationInterface):
    def __init__(self):
        self.frankachain = FrankaChain('panda_link0', ['panda_link1', 'panda_link2', 'panda_link3', 'panda_link4',
                                'panda_link5', 'panda_link6', 'panda_link7'], ['TCP_frame', 'panda_link8', 'panda_hand'], 'panda_arm_model.urdf')
        self.max_vel = pi
        self.n_joints = 7
        self.indJointPos0 = range(0,7)
        self.indTarget =  range(7,10)
        self.indObstaclePos = range(10,13)
        self.indObstacleRadius = 13 
        super().__init__()
        return
    def SetProblemDimensions(self):
        self.nx = self.n_joints 
        self.nu = self.n_joints 
        self.n_stage_params = 1  # dt
        self.n_global_params = 14 # JointPos0, Target, ObstaclePos,  ObstacleRadius
    def Dynamics(self, uk, xk, stage_params, global_params):
        return xk + stage_params[0]*uk
    def StageCost(self, uk, xk, stage_params, global_params):
        ukxk = vertcat(uk, xk)
        H = np.eye(self.nx+self.nu)
        return 0.5*uk.T@uk
    def StageCostFinal(self, xK, stage_params, global_params):
        H = np.eye(self.nx)
        return 0
        # return 0.5 * (xK.T@H@xK)
    def EqConstrInitial(self, uk, xk, stage_params, global_params):
        constr = global_params[self.indJointPos0] - xk 
        self.ngI = constr.shape[0]
        return constr
    def EqConstrFinal(self, xK, stage_params, global_params):
        self.frankachain.set_up_expressions(xK)
        constr = self.frankachain.GetFrameExpression('TCP_frame')[:3,3] - global_params[self.indTarget] 
        self.ngF = constr.shape[0]
        return constr
    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        self.frankachain.set_up_expressions(xk)
        distance_all = self.frankachain.GetAllCapsuleConstraints(global_params[self.indObstaclePos], global_params[self.indObstacleRadius])
        self.distance_all_size = distance_all.shape[0]
        joint_vel = uk
        constr = vertcat(xk, joint_vel, distance_all)
        self.ngIneq = constr.shape[0]
        return constr
    def StageWiseInequalityBounds(self):
        self.lower = np.hstack((self.frankachain.joint_lower, -self.max_vel*np.ones(self.n_joints), np.zeros(self.distance_all_size)))
        self.upper = np.hstack((self.frankachain.joint_upper, self.max_vel*np.ones(self.n_joints), Inf*np.ones(self.distance_all_size)))
        return [self.lower, self.upper]
    def FinalInequality(self, xk, stage_params, global_params):
        constr =  xk
        self.ngIneqF = constr.shape[0]
        return constr
    def FinalInequalityBounds(self):
        return [self.frankachain.joint_lower, self.frankachain.joint_upper]
if __name__ == '__main__':
    robspec = RobSpecCaps()
    K = 10
    dt = 10.0/K
    stage_params = np.matlib.repmat(np.array(dt), 1, K)
    joinstpos0 = np.array([-0.0] +2*[0.0] + [-1.0] + [1.00] + 2*[0.0])
    targetpos = np.array([0.50,0.00,0.50])
    obstablepos = np.array([0.5, 0.00,0.20]) 
    obstaclerad = np.array([0.20])
    global_params =np.hstack((joinstpos0, targetpos, obstablepos, obstaclerad))
    alpha = 0.50
    inits_x = np.matlib.repmat((alpha*robspec.frankachain.joint_lower + (1-alpha)*robspec.frankachain.joint_upper)[:, np.newaxis], 1, K)
    inits_u = np.zeros((robspec.nu, K-1))
    # json + codegen for use in fatrop
    codegen = FatropOCPCodeGenerator(robspec)
    codegen.generate_code('robot.c')
    jsongen = JSONGenerator(robspec)
    jsongen.generate_JSON('robot.json', K, stage_params, global_params, inits_x, inits_u)
    # opti for reference solve
    optibuilder = OptiBuilder(robspec)
    opti = optibuilder.set_up_Opti(K)
    opti.set_value(optibuilder.stage_params_in, stage_params)
    opti.set_value(optibuilder.global_params_in, global_params)
    opti.set_initial(optibuilder.x_vars, inits_x)
    opti.set_initial(optibuilder.u_vars, inits_u)
    # opti.solver('ipopt', {'expand':True}, {'print_level':5, 'hessian_approximation':'limited-memory', 'limited_memory_max_history':1})
    opti.solver('ipopt', {'expand':True}, {'print_level':5})
    res = opti.solve()
    # jsongen.generate_JSON('robot.json', K, stage_params, global_params, np.array(res.value(optibuilder.x_vars)), np.array(res.value(optibuilder.u_vars)))