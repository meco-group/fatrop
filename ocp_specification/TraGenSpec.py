from OCPSpecification import *
import FSDynamics

class TraGenSpec(OCPSpecificationInterface):
    def __init__(self):
        super().__init__()
        self.fsdyns = FSDynamics.FSDynamics()
        self.indp = self.fsdyns.indp
        self.indR = self.fsdyns.indR
        self.ind_params_dt = [0]
        self.ind_params_inv = range(1,4) 
        self.ind_params_end_pos = range(4,7) 
    def SetProblemDimensions(self):
        self.nx = 12
        self.nu = 3
        self.n_stage_params = 1 + 3 + 3 # dt, 3 invariants, 3 end_pos
    def Dynamics(self, uk, xk, stage_params):
        return self.fsdyns.dynamics(uk, xk, stage_params[self.ind_params_dt])
    def StageCost(self, uk, xk, stage_params):
        pass
    def StageCostFinal(self, xK, stage_params):
        pass
    def EqConstrInitial(self, uk, xk, stage_params):
        pass
    def EqConstrFinal(self, xK, stage_params):
        pass
    def StageWiseInequality(self, uk, xk, stage_params):
        pass
    
