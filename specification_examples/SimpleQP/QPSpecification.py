from numpy import Inf
from FatropOCPSpecification import *
from casadi import *

class QPSpecification(OCPSpecificationInterface):
    def __init__(self, w_pos=1, w_rot=1, w_invars=(10**-3)*np.array([1.0, 1.0, 1.0])):
        self.lowerF = np.array([1.0,1.0])
        self.upperF = np.array([2.0,2.0])
        super().__init__()
        return

    def SetProblemDimensions(self):
        self.nx = 2
        self.nu = 2
        self.ngI = 0
        self.ngF =  0 
        self.ngIneq = 0
        self.ngIneqF = 2
        self.n_stage_params = 0  # dt
        self.n_global_params = 0 # endpos 

    def Dynamics(self, uk, xk, stage_params, global_params):
        return xk + uk

    def StageCost(self, uk, xk, stage_params, global_params):
        ukxk = vertcat(uk, xk)
        H = np.eye(self.nx+self.nu)
        return 0.5 * (ukxk.T@H@ukxk)

    def StageCostFinal(self, xK, stage_params, global_params):
        H = np.eye(self.nx)
        return 0.5 * (xK.T@H@xK)
    
    # def EqConstrInitial(self, uk, xk, stage_params, global_params):
    #     pass

    # def EqConstrFinal(self, xK, stage_params, global_params):
    #     pass

    # def StageWiseInequality(self, uk, xk, stage_params, global_params):
    #     pass
    def FinalInequality(self, xk, stage_params, global_params):
        return self.lowerF, xk, self.upperF