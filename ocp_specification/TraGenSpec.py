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
        ## global params
        self.ind_params_R0 = range(9)
        self.ind_params_p0 = range(9,12)
        self.ind_params_RF= range(12,21)
        self.ind_params_pF = range(21,24)
        ## hard_coded params
        self.w_pos = w_pos
        self.w_rot = w_rot
        self.w_invars = w_invars
    def SetProblemDimensions(self):
        self.nx = 12
        self.nu = 3
        self.ngI = 12 
        self.ngF = 12 
        self.ngIneq = 0
        self.n_stage_params = 1 + 3 # dt, 3 invariants
        self.n_global_params = 24 
    def Dynamics(self, uk, xk, stage_params, global_params):
        return self.fsdyns.dynamics(uk, xk, stage_params[self.ind_params_dt])
    def StageCost(self, uk, xk, stage_params, global_params):
        err_invars = self.w_invars*(uk- stage_params[self.ind_params_inv])
        return cas.dot(err_invars, err_invars)
    def StageCostFinal(self, xK, stage_params, global_params):
        return 0
    def EqConstrInitial(self, uk, xk, stage_params, global_params):
        return vertcat(xk[self.indR] - global_params[self.ind_params_R0],xk[self.indp] - global_params[self.ind_params_p0]) 
    def EqConstrFinal(self, xK, stage_params, global_params):
        return vertcat(xK[self.indR] - global_params[self.ind_params_RF],xK[self.indp] - global_params[self.ind_params_pF]) 
    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        return DM.zeros(0), DM.zeros(0), DM.zeros(0) 