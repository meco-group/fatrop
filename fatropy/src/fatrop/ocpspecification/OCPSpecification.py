from abc import abstractmethod
from casadi import *
import json
import numpy as np
import numpy.matlib
from matplotlib.font_manager import json_dump
from typing import List

class BasicOCPInterface:
    def __init__(self):
        # problem dimensions
        self.nx = 0
        self.nu = 0
        self.ngI = 0
        self.ngF = 0
        self.ngIneq = 0
        self.ngIneqF = 0
        self.n_stage_params = 0
        self.n_global_params = 0
        self.SetProblemDimensions()

    @abstractmethod
    def SetProblemDimensions(self):
        pass
    @abstractmethod
    def Dynamics(self, uk, xk, stage_params, global_params):
        pass
    @abstractmethod
    def StageCost(self, uk, xk, stage_params, global_params):
        pass
    @abstractmethod
    def StageCostFinal(self, xK, stage_params, global_params):
        pass
    @abstractmethod
    def EqConstrInitial(self, uk, xk, stage_params, global_params):
        return MX.zeros(0)
    @abstractmethod
    def EqConstrFinal(self, xK, stage_params, global_params):
        return MX.zeros(0)
    @abstractmethod
    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        return MX.zeros(0)
    @abstractmethod
    def StageWiseInequalityBounds(self):
        return np.array([]), np.array([]) 
    @abstractmethod
    def FinalInequality(self, xk, stage_params, global_params):
        return MX.zeros(0)
    @abstractmethod
    def FinalInequalityBounds(self):
        return np.array([]), np.array([]) 
    @abstractmethod
    def DefaultStageParams(self):
        return MX.zeros(self.n_stage_params)
    
class BasicOCPAdapter(FatropOCP):
    def __init__(self, basicocp:BasicOCPInterface, K:int):
        # symbols for generating functions
        uk = MX.sym(basicocp.nu)
        xk = MX.sym(basicocp.nx)
        stage_params = MX.sym(basicocp.n_stage_params)
        global_params = MX.sym(basicocp.n_global_params)
        # initialize functions
        self.AddFunc("Lk", Function([uk, xk, stage_params, global_params], [basicocp.StageCost(uk, xk, stage_params, global_params)]))
        self.AddFunc("LkF", Function([uk, xk, stage_params, global_params], [basicocp.StageCostFinal(xk, stage_params, global_params)]))
        self.AddFunc("dyn", Function([uk, xk, stage_params, global_params], [basicocp.Dynamics(uk, xk, stage_params, global_params)]))
        self.AddFunc("eqI", Function([uk, xk, stage_params, global_params], [basicocp.EqConstrInitial(uk, xk, stage_params, global_params)]))
        self.AddFunc("eqF", Function([uk, xk, stage_params, global_params], [basicocp.EqConstrFinal(xk, stage_params, global_params)]))
        self.AddFunc("ineq", Function([uk, xk, stage_params, global_params], [basicocp.StageWiseInequality(uk, xk, stage_params, global_params)]))
        self.AddFunc("ineqF", Function([uk, xk, stage_params, global_params], [basicocp.FinalInequality(xk, stage_params, global_params)]))
        # initial stage
        self.AddStage(\
        OCPStage(OCPStageDims(basicocp.nu, basicocp.nx, basicocp.n_stage_params, basicocp.ngI, basicocp.ngIneq),\
        OCPStageFuncs("Lk", "dyn", "eqI", "ineq"),\
        OCPStageIneqBounds(basicocp.StageWiseInequalityBounds()[0], basicocp.StageWiseInequalityBounds()[1])))
        # middle stages
        for k in range(1, K-1):
            self.AddStage(\
            OCPStage(OCPStageDims(basicocp.nu, basicocp.nx, basicocp.n_stage_params, basicocp.ngI, basicocp.ngIneq),\
            OCPStageFuncs("Lk", "dyn", "none", "ineq"),\
            OCPStageIneqBounds(basicocp.StageWiseInequalityBounds()[0], basicocp.StageWiseInequalityBounds()[1])))
            # self.AddStage()
        # termminal stage
        self.AddStage(\
        OCPStage(OCPStageDims(basicocp.nu, basicocp.nx, basicocp.n_stage_params, basicocp.ngI, basicocp.ngIneq),\
        OCPStageFuncs("LkF", "none", "eqF", "ineqF"),\
        OCPStageIneqBounds(basicocp.FinalInequalityBounds()[0], basicocp.FinalInequalityBounds()[1])))
        return

class Simulator:
    def __init__(self, ocpspec: BasicOCPInterface):
        self.ocpspec = ocpspec
    def Simulate(self, x0, inputs, global_parms, stage_parms, K):
        nu = self.ocpspec.nu
        nx = self.ocpspec.nx
        n_stage_params = self.ocpspec.n_stage_params
        n_global_params = self.ocpspec.n_global_params
        # make symbols for variables
        u_sym = MX.sym("inputs", nu)
        x_sym = MX.sym("states", nx)
        stage_params_sym = MX.sym("states", n_stage_params)
        global_params_sym = MX.sym("states", n_global_params)
        # make expressions for functions
        dynamics = self.ocpspec.Dynamics(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        # BAbt
        stateskp1 = MX.sym("states_kp1", nx)
        dynamicsF  = Function('dyn',[u_sym, x_sym, stage_params_sym, global_params_sym], [(dynamics)[:]])

        states = np.zeros((nx, K))
        states[:,0] = x0[:]
        for k in range(1, K):
            states[:,[k]] = dynamicsF(inputs[:, k-1], states[:,k-1], stage_parms[:,k-1], global_parms)
        return states


class JSONGenerator:
    def __init__(self, ocpspec: BasicOCPInterface):
        self.ocpspec = ocpspec

    def generate_JSON(self, filename, K, stage_params, global_params, initial_x, initial_u):
        # problem dimensions
        JSONdict = {'nx': self.ocpspec.nx,        'nu': self.ocpspec.nu,        'ngI': self.ocpspec.ngI,        'ngF': self.ocpspec.ngF,
                    'ng_ineq': self.ocpspec.ngIneq, 'ng_ineqF': self.ocpspec.ngIneqF,        'n_stage_params': self.ocpspec.n_stage_params,        'n_global_params': self.ocpspec.n_global_params, 'K': K}
        # stage params
        JSONdict['stage_params'] = stage_params.ravel(order='f').tolist()
        JSONdict['global_params'] = global_params.ravel(order='f').tolist()
        # params
        JSONdict['params'] = stage_params.ravel(order='f').tolist()
        # initial x
        # if initial_x != None:
        JSONdict['initial_x'] = initial_x.ravel(order='f').tolist()

        # else:
        # JSONdict['initial_x'] = np.zeros(self.ocpspec.nx*K).tolist()
        # initial u
        # if initial_u != None:
        JSONdict['initial_u'] = initial_u.ravel(order='f').tolist()
        # else:
        # JSONdict['initial_u'] = np.zeros(self.ocpspec.nu*K).tolist()
        # bounds
        lower, upper = self.ocpspec.StageWiseInequalityBounds()
        lowerF, upperF = self.ocpspec.FinalInequalityBounds()
        JSONdict['lower'] = np.matlib.repmat(lower[:,np.newaxis], 1, K-1).ravel(order='f').tolist()
        JSONdict['upper'] = np.matlib.repmat(upper[:,np.newaxis], 1, K-1).ravel(order='f').tolist()
        JSONdict['lowerF'] =  lowerF.tolist()
        JSONdict['upperF'] =  upperF.tolist()
        print(json_dump(JSONdict, filename))
        return


class FatropOCPCodeGenerator:
    def __init__(self, ocpspec):
        self.ocpspec = ocpspec

    def generate_code(self, filename):
        C = CodeGenerator(filename)
        # get problem dimensions
        nu = self.ocpspec.nu
        nx = self.ocpspec.nx
        n_stage_params = self.ocpspec.n_stage_params
        n_global_params = self.ocpspec.n_global_params
        # make symbols for variables
        u_sym = MX.sym("inputs", nu)
        x_sym = MX.sym("states", nx)
        stage_params_sym = MX.sym("states", n_stage_params)
        global_params_sym = MX.sym("states", n_global_params)
        # make expressions for functions
        dynamics = self.ocpspec.Dynamics(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        Lk = self.ocpspec.StageCost(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        LF = self.ocpspec.StageCostFinal(
            x_sym, stage_params_sym, global_params_sym)
        eqI = self.ocpspec.EqConstrInitial(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        eqF = self.ocpspec.EqConstrFinal(
            x_sym, stage_params_sym, global_params_sym)
        ineq = self.ocpspec.StageWiseInequality(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        ineqF = self.ocpspec.FinalInequality(
            x_sym, stage_params_sym, global_params_sym)
        ngI = eqI.shape[0]
        ngF = eqF.shape[0]
        ngIneq = ineq.shape[0]
        ngIneqF = ineqF.shape[0]
        # make symbols for dual variables
        dual_dyn = MX.sym("d_dyn", nx)
        dual_eqI = MX.sym("d_EqI", ngI)
        dualIneq = MX.sym("dualineq", ngIneq)
        dualIneqF = MX.sym("dualineqF", ngIneqF)
        dual_eqF = MX.sym("d_EqF", ngF)
        obj_scale = MX.sym("obj_scale")
        # BAbt
        stateskp1 = MX.sym("states_kp1", nx)
        BAbt = MX.zeros(nu+nx+1, nx)
        BAbt[:nu+nx,
             :] = jacobian(dynamics, vertcat(u_sym, x_sym)).T
        b = (-stateskp1 + dynamics)[:]
        BAbt[nu+nx, :] = b
        C.add(
            Function("BAbt", [stateskp1, u_sym, x_sym, stage_params_sym, global_params_sym], [densify(BAbt)]).expand())
        # b
        C.add(Function("bk", [stateskp1, u_sym,
                              x_sym, stage_params_sym, global_params_sym], [densify(b)]).expand())
        # RSQrqtI
        RSQrqtI = MX.zeros(nu+nx+1, nu + nx)
        [RSQI, rqI] = hessian(Lk, vertcat(u_sym, x_sym))
        RSQIGN = RSQI
        rqlagI = rqI
        if ngI > 0:
            [H, h]= hessian(dual_eqI.T@eqI, vertcat(u_sym, x_sym))
            RSQI += H
            rqlagI += h
        [H,h] = hessian(dual_dyn.T@dynamics,
                        vertcat(u_sym, x_sym))
        RSQI += H
        rqlagI += h
        
        if ngIneq > 0:
            [H,h] = hessian(dualIneq.T@ineq,
                            vertcat(u_sym, x_sym))
            RSQI += H
            rqlagI += h
        RSQrqtI[:nu+nx, :] = RSQI
        RSQrqtI[nu+nx, :] = rqlagI[:]
        C.add(Function("RSQrqtI", [obj_scale, u_sym,
              x_sym, dual_dyn, dual_eqI, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqtI)]).expand())
        RSQrqtI[:nu+nx, :] = RSQIGN
        RSQrqtI[nu+nx, :] = rqlagI[:]
        C.add(Function("RSQrqtIGN", [obj_scale, u_sym,
              x_sym, dual_dyn, dual_eqI, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqtI)]).expand())
        # rqI
        C.add(Function("rqI", [obj_scale,
              u_sym, x_sym, stage_params_sym, global_params_sym], [densify(rqI)]).expand())
        # RSQrqt
        RSQrqt = MX.zeros(nu+nx+1, nu + nx)
        [RSQ, rq] = hessian(Lk, vertcat(u_sym, x_sym))
        RSQGN = RSQ
        rqlag = rq
        [H,h]= hessian(dual_dyn.T@dynamics,
                       vertcat(u_sym, x_sym))
        RSQ += H
        rqlag +=h

        if ngIneq > 0:
            [H,h] = hessian(dualIneq.T@ineq,
                           vertcat(u_sym, x_sym))
            RSQ += H
            rqlag +=h
        RSQrqt[:nu+nx, :] = RSQ
        RSQrqt[nu+nx, :] = rqlag[:]
        C.add(Function("RSQrqt", [obj_scale, u_sym, x_sym,
              dual_dyn, dual_eqI, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqt)]).expand())
        RSQrqt[:nu+nx, :] = RSQGN
        RSQrqt[nu+nx, :] = rqlag[:]
        C.add(Function("RSQrqtGN", [obj_scale, u_sym, x_sym,
              dual_dyn, dual_eqI, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqt)]).expand())
        # rqF
        C.add(Function("rqk", [obj_scale,
              u_sym, x_sym, stage_params_sym, global_params_sym], [densify(rq)]).expand())
        # Lk
        C.add(Function("Lk", [obj_scale, u_sym,
              x_sym, stage_params_sym, global_params_sym], [densify(Lk)]).expand())
        # RSQrqtF
        RSQrqtF = MX.zeros(nx+1, nx)
        [RSQF, rqF] = hessian(LF, vertcat(x_sym))
        RSQFGN = RSQF
        rqlagF = rqF
        if ngF > 0:
            [H, h]= hessian(dual_eqF.T@eqF,
                            vertcat(x_sym))
            RSQF += H
            rqlagF += h
        if ngIneqF > 0:
            [H,h] = hessian(dualIneqF.T@ineqF,
                           vertcat(x_sym))
            RSQF += H
            rqlagF += h
        # if ngIneq>-1:
        #     RSQF += hessian(dualIneq.T@ineq, vertcat(u_sym, x_sym))[-1]
        RSQrqtF[:nx, :] = RSQF
        RSQrqtF[nx, :] = rqlagF[:]
        C.add(Function("RSQrqtF", [obj_scale, u_sym, x_sym,
              dual_dyn, dual_eqF, dualIneqF, stage_params_sym, global_params_sym], [densify(RSQrqtF)]).expand())
        RSQrqtF[:nx, :] = RSQFGN
        RSQrqtF[nx, :] = rqlagF[:]
        C.add(Function("RSQrqtFGN", [obj_scale, u_sym, x_sym,
              dual_dyn, dual_eqF, dualIneqF, stage_params_sym, global_params_sym], [densify(RSQrqtF)]).expand())
        # rqF
        C.add(Function("rqF", [obj_scale,
              u_sym, x_sym, stage_params_sym, global_params_sym], [densify(rqF)]).expand())
        # LF
        C.add(Function("LF", [obj_scale, u_sym,
              x_sym, stage_params_sym, global_params_sym], [densify(LF)]).expand())
        # GgtI
        GgtI = MX.zeros(nu+nx+1, ngI)
        GgtI[:nu+nx,
             :] = jacobian(eqI, vertcat(u_sym, x_sym)).T
        GgtI[nu+nx, :] = eqI[:].T
        C.add(Function(
            "GgtI", [u_sym, x_sym, stage_params_sym, global_params_sym], [densify(GgtI)]).expand())
        # g_I
        C.add(Function("gI", [u_sym, x_sym, stage_params_sym,
              global_params_sym], [densify(eqI[:])]).expand())
        # GgtF
        GgtF = MX.zeros(nx+1, ngF)
        GgtF[:nx, :] = jacobian(eqF, x_sym).T
        GgtF[nx, :] = eqF[:].T
        C.add(Function(
            "GgtF", [u_sym, x_sym, stage_params_sym, global_params_sym], [densify(GgtF)]).expand())
        # g_F
        C.add(Function("gF", [u_sym, x_sym, stage_params_sym,
              global_params_sym], [densify(eqF[:])]).expand())
        # Ggineqt
        Ggineqt = MX.zeros(nu+nx+1, ngIneq)
        Ggineqt[:nu+nx,
                :] = jacobian(ineq, vertcat(u_sym, x_sym)).T
        Ggineqt[nu+nx, :] = ineq[:].T
        C.add(Function("Ggineqt", [u_sym,
              x_sym, stage_params_sym, global_params_sym], [densify(Ggineqt)]).expand())
        C.add(Function("gineq", [u_sym, x_sym, stage_params_sym, global_params_sym], [
              densify(ineq[:])]).expand())
        # GgineqFt
        GgineqFt = MX.zeros(nx+1, ngIneqF)
        GgineqFt[:nx,
                :] = jacobian(ineqF, vertcat(x_sym)).T
        GgineqFt[nx, :] = ineqF[:].T
        C.add(Function("GgineqFt", [
              x_sym, stage_params_sym, global_params_sym], [densify(GgineqFt)]).expand())
        C.add(Function("gineqF", [x_sym, stage_params_sym, global_params_sym], [
              densify(ineqF[:])]).expand())
        C.add(Function("default_stage_params", [], [
              densify(self.ocpspec.DefaultStageParams())]).expand())

        C.generate()
        return

class OptiBuilder:
    def __init__(self, ocpspec: BasicOCPInterface):
        self.ocpspec = ocpspec

    def set_up_Opti(self, K: int):
        # all scales are set to 1.0
        # get problem dimensions
        nu = self.ocpspec.nu
        nx = self.ocpspec.nx
        self.opti = Opti()
        self.stage_params_in = self.opti.parameter(
            self.ocpspec.n_stage_params, K)
        self.global_params_in = self.opti.parameter(
            self.ocpspec.n_global_params, 1)
        self.N_vars = K*nx + (K-1)*nu
        self.x_vars = self.opti.variable(nx, K)
        self.u_vars = self.opti.variable(nu, K-1)
        # make symbols for variables
        u_sym = MX.sym("inputs", nu)
        x_sym = MX.sym("states", nx)
        stage_params_sym = MX.sym("stageparams", self.ocpspec.n_stage_params)
        global_params_sym = MX.sym("stageparams", self.ocpspec.n_global_params)
        obj_scale = MX.sym("obj_scale", 1)
        # make expressions for functions
        dynamics = self.ocpspec.Dynamics(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        Lk = self.ocpspec.StageCost(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        LF = self.ocpspec.StageCostFinal(
            x_sym, stage_params_sym, global_params_sym)
        eqI = self.ocpspec.EqConstrInitial(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        eqF = self.ocpspec.EqConstrFinal(
            x_sym, stage_params_sym, global_params_sym)
        ineq = self.ocpspec.StageWiseInequality(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        lower, upper = self.ocpspec.StageWiseInequalityBounds()
        ineqF = self.ocpspec.FinalInequality(
            x_sym, stage_params_sym, global_params_sym)
        lowerF, upperF = self.ocpspec.FinalInequalityBounds()
        ngI = eqI.shape[0]
        ngF = eqF.shape[0]
        ngIneq = ineq.shape[0]
        ngIneqF = ineqF.shape[0]
        Lkf = Function(
            "Lk", [obj_scale, x_sym, u_sym, stage_params_sym, global_params_sym], [Lk])
        LkFf = Function(
            "LF", [obj_scale, x_sym, stage_params_sym, global_params_sym], [LF])
        stateskp1 = MX.sym("stateskp1", nx)
        Dynamcisf = Function(
            "F", [stateskp1, x_sym, u_sym, stage_params_sym, global_params_sym], [-stateskp1 + dynamics])
        if ngI > 0:
            EqIf = Function(
                "eqI", [x_sym, u_sym, stage_params_sym, global_params_sym], [eqI])
            self.opti.subject_to(
                EqIf(self.x_vars[:, 0], self.u_vars[:, 0], self.stage_params_in[:, 0], self.global_params_in) == 0.0)
        if ngF > 0:
            EqFf = Function(
                "eqF", [x_sym, stage_params_sym, global_params_sym], [eqF])
            self.opti.subject_to(EqFf(
                self.x_vars[:, K-1], self.stage_params_in[:, K-1], self.global_params_in) == 0.0)
        J = 0
        for k in range(K-1):
            J += Lkf(1.0, self.x_vars[:, k], self.u_vars[:, k],
                     self.stage_params_in[:, k], self.global_params_in)
            self.opti.subject_to(
                Dynamcisf(self.x_vars[:, k+1], self.x_vars[:, k], self.u_vars[:, k], self.stage_params_in[:, k], self.global_params_in) == 0.0)
        if ngIneq > 0:
            Ineqf = Function(
                "ineqf", [u_sym, x_sym, stage_params_sym, global_params_sym], [ineq])
            for k in range(K-1):
                for i in range(ngIneq):
                    if lower[i] == -inf:
                        # self.opti.subject_to(Ineqf(self.u_sym, self.x_sym)[i]< self.upper[i])
                        self.opti.subject_to(upper[i] > Ineqf(
                            self.u_vars[:, k], self.x_vars[:, k], self.stage_params_in[:, k], self.global_params_in)[i])
                    elif upper[i] == inf:
                        self.opti.subject_to(lower[i] < Ineqf(
                            self.u_vars[:, k], self.x_vars[:, k], self.stage_params_in[:, k], self.global_params_in)[i])
                    else:
                        self.opti.subject_to(lower[i] < (Ineqf(
                            self.u_vars[:,k], self.x_vars[:,k], self.stage_params_in[:, k], self.global_params_in)[i] < upper[i]))
        if ngIneqF > 0:
            IneqFf = Function(
                "ineqFf", [x_sym, stage_params_sym, global_params_sym], [ineqF])
            for i in range(ngIneqF):
                if lowerF[i] == -inf:
                    # self.opti.subject_to(IneqFf(self.u_sym, self.x_sym)[i]< self.upper[i])
                    self.opti.subject_to(upperF[i] > IneqFf(
                       self.x_vars[:, K-1], self.stage_params_in[:, K-1], self.global_params_in)[i])
                elif upperF[i] == inf:
                    self.opti.subject_to(lowerF[i] < IneqFf(
                        self.x_vars[:, K-1], self.stage_params_in[:, K-1], self.global_params_in)[i])
                else:
                    self.opti.subject_to(lowerF[i] < (IneqFf(
                        self.x_vars[:,K-1], self.stage_params_in[:, K-1], self.global_params_in)[i] < upperF[i]))
        J += LkFf(1.0, self.x_vars[:, K-1],
                  self.stage_params_in[:, K-1], self.global_params_in)
        self.opti.minimize(J)
        return self.opti
    pass


# # This class can be used to specify OCP's with intial and final constraints and generate necessary code for it.
# class OptimalControlProblem:
#     def __init__(self):
#         self.ngI = 0
#         self.ngF = 0
#         self.ngIneq = 0
#         self.eqI = MX.sym("EqI", 0)
#         self.dual_eqI = MX.sym("d_EqI", 0)
#         self.dualIneq = MX.sym("dualineq", 0)
#         self.eqF = MX.sym("d_EqI", 0)
#         self.dual_eqF = MX.sym("d_EqF", 0)
#         self.obj_scale = MX.sym("obj_scale")
#         self.ineq = MX.sym("d_EqI", 0)

#     def get_states(self, nx):
#         self.nx = nx
#         self.x_sym = MX.sym('states', nx)
#         return self.x_sym

#     def get_inputs(self, nu):
#         self.nu = nu
#         self.u_sym = MX.sym('inputs', nu)
#         return self.u_sym

#     def set_dynamics(self, xkp1):
#         self.dual_dyn = MX.sym("dual_dyn", self.nx)
#         self.dynamics = xkp1

#     def set_stagecost(self, Lk):
#         self.Lk = self.obj_scale*Lk

#     def set_stagecostFinal(self, LF):
#         self.LF = self.obj_scale*LF

#     def set_eq_initial(self, eq):
#         self.ngI = eq.shape[0]
#         self.dual_eqI = MX.sym("dual_eqI", self.ngI)
#         self.eqI = eq

#     def set_eq_final(self, eq):
#         self.ngF = eq.shape[0]
#         self.dual_eqF = MX.sym("dual_eqF", self.ngF)
#         self.eqF = eq

#     def set_ineq(self, ineq, lower, upper):
#         # todo no ineq on final stage x
#         self.ngIneq = ineq.shape[0]
#         self.dualIneq = MX.sym("dual_Ineq", self.ngIneq)
#         self.ineq = ineq
#         self.lower = lower
#         self.upper = upper

#     def set_up_Opti(self, K):
#         # all scales are set to 1.0
#         nu = self.nu
#         nx = self.nx
#         self.opti = Opti()
#         self.N_vars = K*nx + (K-1)*nu
#         self.opti_vars = self.opti.variable(self.N_vars)
#         self.x_vars = MX.zeros(nx, K)
#         self.u_vars = MX.zeros(nu, K)
#         for k in range(K-1):
#             offs = k*(nu+nx)
#             self.u_vars[:, k] = self.opti_vars[offs:offs+nu]
#             self.x_vars[:, k] = self.opti_vars[offs+nu:offs+nu+nx]
#         offs = (K-1)*(nu+nx)
#         self.x_vars[:, K-1] = self.opti_vars[offs:offs+nx]
#         Lkf = Function(
#             "Lk", [self.obj_scale, self.x_sym, self.u_sym], [self.Lk])
#         LkFf = Function("LF", [self.obj_scale, self.x_sym], [self.LF])
#         stateskp1 = MX.sym("stateskp1", nx)
#         Dynamcisf = Function(
#             "F", [stateskp1, self.x_sym, self.u_sym], [-stateskp1 + self.dynamics])
#         if self.ngI > 0:
#             EqIf = Function("eqI", [self.x_sym, self.u_sym], [self.eqI])
#             self.opti.subject_to(
#                 EqIf(self.x_vars[:, 0], self.u_vars[:, 0]) == 0.0)
#         if self.ngF > 0:
#             EqFf = Function("eqF", [self.x_sym], [self.eqF])
#             self.opti.subject_to(EqFf(self.x_vars[:, K-1]) == 0.0)
#         J = 0
#         for k in range(K-1):
#             J += Lkf(1.0, self.x_vars[:, k], self.u_vars[:, k])
#             self.opti.subject_to(
#                 Dynamcisf(self.x_vars[:, k+1], self.x_vars[:, k], self.u_vars[:, k]) == 0.0)
#         if self.ngIneq > 0:
#             Ineqf = Function("ineqf", [self.u_sym, self.x_sym], [self.ineq])
#             for k in range(K-1):
#                 for i in range(self.ngIneq):
#                     if self.lower[i] == -inf:
#                         # self.opti.subject_to(Ineqf(self.u_sym, self.x_sym)[i]< self.upper[i])
#                         self.opti.subject_to(self.upper[i] > Ineqf(
#                             self.u_vars[:, k], self.x_vars[:, k])[i])
#                     elif self.upper[i] == inf:
#                         self.opti.subject_to(self.lower[i] < Ineqf(
#                             self.u_vars[:, k], self.x_vars[:, k])[i])
#                     else:
#                         self.opti.subject_to(self.lower[i] < Ineqf(
#                             self.u_sym, self.x_sym)[i] < self.upper[i])
#         J += LkFf(1.0, self.x_vars[:, K-1])
#         self.opti.minimize(J)
#         return self.opti

#     def generate_code(self, filename):
#         C = CodeGenerator(filename)
#         # BAbt
#         stateskp1 = MX.sym("states_kp1", self.nx)
#         BAbt = MX.zeros(self.nu+self.nx+1, self.nx)
#         BAbt[:self.nu+self.nx,
#              :] = jacobian(self.dynamics, vertcat(self.u_sym, self.x_sym)).T
#         b = (-stateskp1 + self.dynamics)[:]
#         BAbt[self.nu+self.nx, :] = b
#         C.add(
#             Function("BAbt", [stateskp1, self.u_sym, self.x_sym], [densify(BAbt)]))
#         # b
#         C.add(Function("bk", [stateskp1, self.u_sym,
#                               self.x_sym], [densify(b)]))
#         # RSQrqtI
#         RSQrqtI = MX.zeros(self.nu+self.nx+1, self.nu + self.nx)
#         [RSQI, rqI] = hessian(self.Lk, vertcat(self.u_sym, self.x_sym))
#         if self.ngI > 0:
#             RSQI += hessian(self.dual_eqI.T@self.eqI,
#                             vertcat(self.u_sym, self.x_sym))[0]
#         RSQI += hessian(self.dual_dyn.T@self.dynamics,
#                         vertcat(self.u_sym, self.x_sym))[0]
#         if self.ngIneq > 0:
#             RSQI += hessian(self.dualIneq.T@self.ineq,
#                             vertcat(self.u_sym, self.x_sym))[0]
#         RSQrqtI[:self.nu+self.nx, :] = RSQI
#         RSQrqtI[self.nu+self.nx, :] = rqI[:]
#         C.add(Function("RSQrqtI", [self.obj_scale, self.u_sym,
#               self.x_sym, self.dual_dyn, self.dual_eqI, self.dualIneq], [densify(RSQrqtI)]))
#         rqI
#         C.add(Function("rqI", [self.obj_scale,
#               self.u_sym, self.x_sym], [densify(rqI)]))
#         # RSQrqt
#         RSQrqt = MX.zeros(self.nu+self.nx+1, self.nu + self.nx)
#         [RSQ, rq] = hessian(self.Lk, vertcat(self.u_sym, self.x_sym))
#         RSQ += hessian(self.dual_dyn.T@self.dynamics,
#                        vertcat(self.u_sym, self.x_sym))[0]
#         if self.ngIneq > 0:
#             RSQ += hessian(self.dualIneq.T@self.ineq,
#                            vertcat(self.u_sym, self.x_sym))[0]
#         RSQrqt[:self.nu+self.nx, :] = RSQ
#         RSQrqt[self.nu+self.nx, :] = rq[:]
#         C.add(Function("RSQrqt", [self.obj_scale, self.u_sym, self.x_sym,
#               self.dual_dyn, self.dual_eqI, self.dualIneq], [densify(RSQrqt)]))
#         # rqF
#         C.add(Function("rqk", [self.obj_scale,
#               self.u_sym, self.x_sym], [densify(rq)]))
#         # Lk
#         C.add(Function("Lk", [self.obj_scale, self.u_sym,
#               self.x_sym], [densify(self.Lk)]))
#         # RSQrqtF
#         RSQrqtF = MX.zeros(self.nx+1, self.nx)
#         [RSQF, rqF] = hessian(self.LF, vertcat(self.x_sym))
#         if self.ngF > 0:
#             RSQF += hessian(self.dual_eqF.T@self.eqF,
#                             vertcat(self.x_sym))[0]
#         # if self.ngIneq>-1:
#         #     RSQF += hessian(self.dualIneq.T@self.ineq, vertcat(self.u_sym, self.x_sym))[-1]
#         RSQrqtF[:self.nx, :] = RSQF
#         RSQrqtF[self.nx, :] = rqF[:]
#         C.add(Function("RSQrqtF", [self.obj_scale, self.u_sym, self.x_sym,
#               self.dual_dyn, self.dual_eqF, self.dualIneq], [densify(RSQrqtF)]))
#         # rqF
#         C.add(Function("rqF", [self.obj_scale,
#               self.u_sym, self.x_sym], [densify(rqF)]))
#         # LF
#         C.add(Function("LF", [self.obj_scale, self.u_sym,
#               self.x_sym], [densify(self.LF)]))
#         # GgtI
#         GgtI = MX.zeros(self.nu+self.nx+1, self.ngI)
#         GgtI[:self.nu+self.nx,
#              :] = jacobian(self.eqI, vertcat(self.u_sym, self.x_sym)).T
#         GgtI[self.nu+self.nx, :] = self.eqI[:].T
#         C.add(Function("GgtI", [self.u_sym, self.x_sym], [densify(GgtI)]))
#         # g_I
#         C.add(Function("gI", [self.u_sym, self.x_sym], [densify(self.eqI[:])]))
#         # GgtF
#         GgtF = MX.zeros(self.nx+1, self.ngF)
#         GgtF[:self.nx, :] = jacobian(self.eqF, self.x_sym).T
#         GgtF[self.nx, :] = self.eqF[:].T
#         C.add(Function("GgtF", [self.u_sym, self.x_sym], [densify(GgtF)]))
#         # g_F
#         C.add(Function("gF", [self.u_sym, self.x_sym], [densify(self.eqF[:])]))
#         # Ggineqt
#         Ggineqt = MX.zeros(self.nu+self.nx+1, self.ngIneq)
#         Ggineqt[:self.nu+self.nx,
#                 :] = jacobian(self.ineq, vertcat(self.u_sym, self.x_sym)).T
#         Ggineqt[self.nu+self.nx, :] = self.ineq[:].T
#         C.add(Function("Ggineqt", [self.u_sym,
#               self.x_sym], [densify(Ggineqt)]))
#         C.add(Function("gineq", [self.u_sym, self.x_sym], [
#               densify(self.ineq[:])]))
#         C.generate()
#         return
class OptiCodeGenerator:
    def __init__(self, opti_in:Opti, filename):
        self.opti_in = opti_in
        self.filename = filename
        self.ng = opti_in.ng
        self.np = opti_in.np
        self.nx = opti_in.nx
    def GenerateCode(self):
        lag = self.opti_in.j + self.opti_in.lam_g.T@self.opti_in.g
        constr = self.opti_in.g
        hess_lag = hessian(lag, self.opti_in.x)[0]
        grad_f = jacobian(self.opti_in.j, self.opti_in.x)
        jac_g = jacobian(constr, self.opti_in.x)
        cg = CodeGenerator(self.filename + '.c')
        cg.add(Function('hess_lag', [self.opti_in.x, self.opti_in.p], [hess_lag]).expand())
        cg.add(Function('grad_f', [self.opti_in.x, self.opti_in.p], [grad_f]).expand())
        cg.add(Function('jac_g', [self.opti_in.x, self.opti_in.p], [jac_g]).expand())
        cg.generate()
    def GetOpti():
        pass

##### NOTE this code currently only supports specification of OCP's of the 'basic' type. 
##### fatrop is able to solve a more general class of ocp's the code from here is ongoing
##### work on supporting specification of very general optimal control problems.

class OCPStageDims:
    def __init__(self, nu, nx, nstage_params, ng, ngineq):
        self.nu = nu
        self.nx = nx
        self.nstage_params = nstage_params
        self.ng= ng
        self.ngineq = ngineq
class OCPStageIneqBounds:
    def __init__(self, lower_bounds, upper_bounds):
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds

class OCPStageFuncs:
    def __init__(self, costf:str, dynamicsf:str, eqcf:str, ineqcf:str):
        self.costf = costf
        self.dynamicsf = dynamicsf
        self.eqcf = eqcf
        self.ineqf = ineqcf

class OCPStage:
    def __init__(self, ocpstagedims:OCPStageDims, ocpstagefuncs:OCPStageFuncs, ocpstageineqbounds:OCPStageIneqBounds):
        self.ocpstagedims = ocpstagedims
        self.ocpstagefuncs = ocpstagefuncs
        self.stageineqbounds = ocpstageineqbounds
    
class FatropOCP:
    def __init__(self):
        self.stages:List[OCPStage] = []
        self.functionDict = {}
    def AddStage(self, ocpstage:OCPStage):
        self.stages.append(ocpstage)
    def AddFunc(self, name:str, func:Function):
        self.functionDict[name] = func

class FatropFunctionGenerator:
    ## function name identifiers format
    ## Lk -> 001_
    ## Dynamics -> 002_
    ## eq -> 003_
    ## ineq -> 004_
    ## BAbt -> 005_ + dyn
    ## b -> 006_ + dyn
    ## RSQrq -> 007_ + Lk + dyn + eqs + ineqs
    ## rq -> 008_ + Lk
    ## RSQrq_GN -> 009_ + Lk + dyn + eqs + ineqs
    ## Ggt -> 010_ + eq
    ## g -> 011_ + eq
    ## Ggineqt -> 012_ + ineq
    ## gineq -> 013_ + ineq
    def __init__(self, ocp:FatropOCP):
        self.functions = {}
        K = len(ocp.stages)
        pass
    def GenerateStage(self, stage:OCPStage, functionsdict):
        nx = stage.ocpstagedims.nx
        nu = stage.ocpstagedims.nu
        nstage_params = stage.ocpstagedims.nstage_params
        ng = stage.ocpstagedims.ng
        ngineq = stage.ocpstagedims.ngineq
        ## Lk -> 001_
        name = "001_" + stage.ocpstagefuncs.costf
        if name not in self.functions:
            self.functions[name] = functionsdict[stage.ocpstagefuncs.costf]
        ## Dynamics -> 002_
        ## eq -> 003_
        ## ineq -> 004_
        ## BAbt -> 005_ + dyn
        ## b -> 006_ + dyn
        ## RSQrq -> 007_ + Lk + dyn + eqs + ineqs
        ## rq -> 008_ + Lk
        ## RSQrq_GN -> 009_ + Lk + dyn + eqs + ineqs
        ## Ggt -> 010_ + eq
        ## g -> 011_ + eq
        ## Ggineqt -> 012_ + ineq
        ## gineq -> 013_ + ineq
        pass
    def Generate(self):
        pass