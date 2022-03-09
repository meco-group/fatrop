from numpy import Inf
from FatropOCPSpecification import *
from casadi import *


class QPSpecification(OCPSpecificationInterface):
    def __init__(self, w_pos=1, w_rot=1, w_invars=(10**-3)*np.array([1.0, 1.0, 1.0])):
        self.lowerF = np.array([1.0])
        self.upperF = np.array([inf])
        super().__init__()
        return

    def SetProblemDimensions(self):
        self.nx = 2
        self.nu = 2
        self.ngI = 0
        self.ngF = 1
        self.ngIneq = 1
        self.ngIneqF = 1
        self.n_stage_params = 0  # dt
        self.n_global_params = 0  # endpos

    def Dynamics(self, uk, xk, stage_params, global_params):
        return xk + uk

    def StageCost(self, uk, xk, stage_params, global_params):
        ukxk = vertcat(uk, xk)
        H = np.eye(self.nx+self.nu)
        return 0.5 * (ukxk.T@H@ukxk) + 10*sum1(xk) + sum1(uk)
    def EqConstrFinal(self, xK, stage_params, global_params):
        return xK[0] 

    def StageCostFinal(self, xK, stage_params, global_params):
        H = np.eye(self.nx)
        return 0.5 * (xK.T@H@xK) + 10*sum1(xK)
    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        return xk[1]
    def StageWiseInequalityBounds(self):
        return np.array([10.0]), np.array([inf]) 

    def FinalInequality(self, xk, stage_params, global_params):
        return xk[1]

    def FinalInequalityBounds(self):
        return self.lowerF, self.upperF
