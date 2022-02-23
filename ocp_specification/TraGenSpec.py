from OCPSpecification import *
import FSDynamics
import casadi as cas

class TraGenSpec(OCPSpecificationInterface):
    def __init__(self, w_pos = 1, w_rot = 1, w_invars = (10**-3)*np.array([1.0, 1.0, 1.0])):
        super().__init__()
        self.fsdyns = FSDynamics.FSDynamics()
        self.indp = self.fsdyns.indp
        self.indR = self.fsdyns.indR
        self.ind_params_dt = [0]
        self.ind_params_inv = range(1,4) 
        self.ind_params_end_pos = range(4,7) 
        self.w_pos = w_pos
        self.w_rot = w_rot
        self.w_invars = w_invars
    def SetProblemDimensions(self):
        self.nx = 12
        self.nu = 3
        self.n_stage_params = 1 + 3 + 3 # dt, 3 invariants, 3 end_pos
    def Dynamics(self, uk, xk, stage_params):
        return self.fsdyns.dynamics(uk, xk, stage_params[self.ind_params_dt])
    def StageCost(self, uk, xk, stage_params):
        err_invars = self.w_invars*(uk- stage_params[self.ind_params_inv])
        return cas.dot(err_invars, err_invars)
    def StageCostFinal(self, xK, stage_params):
        return 0
    def EqConstrInitial(self, uk, xk, stage_params):
        pass
    def EqConstrFinal(self, xK, stage_params):
        pass
    def StageWiseInequality(self, uk, xk, stage_params):
        pass
    
