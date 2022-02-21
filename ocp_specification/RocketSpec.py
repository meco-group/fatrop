from OCPSpecification import *
import RocketDynamics

class RocketSpec(OCPSpecificationInterface):
    def __init__(self):
       self.ind_pos = list(range(0,2))
       self.ind_rot = list(range(2,6))
       self.ind_vel = list(range(6,8))
       self.ind_omega = [8]
       super().__init__()
    def SetProblemDimensions(self):
        self.nx = 9
        self.nu = 2
        self.rocketdynamics = RocketDynamics.RockDyns()
    def Dynamics(self, uk, xk):
        return self.rocketdynamics.dynamics(uk, xk, 0.1)
    def StageCost(self, uk, xk):
        return xk.T@xk + sum1(xk) + uk.T@uk
    def StageCostFinal(self, xK):
        return xK.T@xK
    def EqConstrInitial(self, uk, xk):
        return vertcat(xk[self.ind_pos],xk[self.ind_vel])
    def EqConstrFinal(self, xK):
        return vertcat(xK[self.ind_pos]-10*DM.ones(2,1),xK[self.ind_vel])
    def StageWiseInequality(self, uk, xk):
        return [0], uk[0], [inf]