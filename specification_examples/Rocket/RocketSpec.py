from FatropOCPSpecification import *

import RocketDynamics

class RocketSpec(OCPSpecificationInterface):
    def __init__(self):
       super().__init__()
       self.ind_pos = list(range(0,2))
       self.ind_rot = list(range(2,6))
       self.ind_vel = list(range(6,8))
       self.ind_omega = [8]
    def SetProblemDimensions(self):
        self.nx = 9
        self.nu = 2
        self.ngIneq = 1
        self.ngIneqF = 2
        self.ngF = 2
        self.ngI = 4
        self.rocketdynamics = RocketDynamics.RockDyns()
    def Dynamics(self, uk, xk, stage_params, global_params):
        return self.rocketdynamics.dynamics(uk, xk, 0.1)
    def StageCost(self, uk, xk, stage_params, global_params):
        return xk.T@xk + sum1(xk) + uk.T@uk
    def StageCostFinal(self, xK, stage_params, global_params):
        return xK.T@xK
    def EqConstrInitial(self, uk, xk, stage_params, global_params):
        return vertcat(xk[self.ind_pos],xk[self.ind_vel])
    def EqConstrFinal(self, xK, stage_params, global_params):
        return vertcat(xK[self.ind_pos]-10*DM.ones(2,1))
    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        return uk[0]
    def StageWiseInequalityBounds(self):
        return [np.array([0.0]), np.array([inf])]
    def FinalInequality(self, xk, stage_params, global_params):
        return xk[self.ind_vel]
    def FinalInequalityBounds(self):
        return [np.array([0, 0]), np.array([1e-1, 1e-1])]
