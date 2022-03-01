from abc import abstractmethod
from casadi import *
import json

from matplotlib.font_manager import json_dump


class OCPSpecificationInterface:
    def __init__(self):
        # problem dimensions
        self.nx = 0
        self.nu = 0
        self.ngI = 0
        self.ngF = 0
        self.ngIneq = 0
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
        return SX.zeros(0)

    @abstractmethod
    def EqConstrFinal(self, xK, stage_params):
        return SX.zeros(0)

    @abstractmethod
    def StageWiseInequality(self, uk, xk, stage_params, global_params):
        return SX.zeros(0)

    @abstractmethod
    def DefaultStageParams(self):
        return SX.zeros(self.n_stage_params)


class JSONGenerator:
    def __init__(self, ocpspec):
        self.ocpspec = ocpspec

    def generate_JSON(self, filename, K, stage_params, global_params, initial_x, initial_u, lower =np.array([]), upper = np.array([])):
        # problem dimensions
        JSONdict = {'nx': self.ocpspec.nx,        'nu': self.ocpspec.nu,        'ngI': self.ocpspec.ngI,        'ngF': self.ocpspec.ngF,
                    'ng_ineq': self.ocpspec.ngIneq,        'n_stage_params': self.ocpspec.n_stage_params,        'n_global_params': self.ocpspec.n_global_params, 'K': K}
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
        JSONdict['lower'] = lower.ravel(order='f').tolist()
        JSONdict['upper'] = upper.ravel(order='f').tolist()
        print(json_dump(JSONdict, "test.json"))
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
        u_sym = SX.sym("inputs", nu)
        x_sym = SX.sym("states", nx)
        stage_params_sym = SX.sym("states", n_stage_params)
        global_params_sym = SX.sym("states", n_global_params)
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
            u_sym, x_sym, stage_params_sym, global_params_sym)[1]
        ngI = eqI.shape[0]
        ngF = eqF.shape[0]
        ngIneq = ineq.shape[0]
        # make symbols for dual variables
        dual_dyn = SX.sym("d_dyn", nx)
        dual_eqI = SX.sym("d_EqI", ngI)
        dualIneq = SX.sym("dualineq", ngIneq)
        dual_eqF = SX.sym("d_EqF", ngF)
        obj_scale = SX.sym("obj_scale")
        # BAbt
        stateskp1 = SX.sym("states_kp1", nx)
        BAbt = SX.zeros(nu+nx+1, nx)
        BAbt[:nu+nx,
             :] = jacobian(dynamics, vertcat(u_sym, x_sym)).T
        b = (-stateskp1 + dynamics)[:]
        BAbt[nu+nx, :] = b
        C.add(
            Function("BAbt", [stateskp1, u_sym, x_sym, stage_params_sym, global_params_sym], [densify(BAbt)]))
        # b
        C.add(Function("bk", [stateskp1, u_sym,
                              x_sym, stage_params_sym, global_params_sym], [densify(b)]))
        # RSQrqtI
        RSQrqtI = SX.zeros(nu+nx+1, nu + nx)
        [RSQI, rqI] = hessian(Lk, vertcat(u_sym, x_sym))
        if ngI > 0:
            RSQI += hessian(dual_eqI.T@eqI, vertcat(u_sym, x_sym))[0]
        RSQI += hessian(dual_dyn.T@dynamics,
                        vertcat(u_sym, x_sym))[0]
        if ngIneq > 0:
            RSQI += hessian(dualIneq.T@ineq,
                            vertcat(u_sym, x_sym))[0]
        RSQrqtI[:nu+nx, :] = RSQI
        RSQrqtI[nu+nx, :] = rqI[:]
        C.add(Function("RSQrqtI", [obj_scale, u_sym,
              x_sym, dual_dyn, dual_eqI, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqtI)]))
        rqI
        C.add(Function("rqI", [obj_scale,
              u_sym, x_sym, stage_params_sym, global_params_sym], [densify(rqI)]))
        # RSQrqt
        RSQrqt = SX.zeros(nu+nx+1, nu + nx)
        [RSQ, rq] = hessian(Lk, vertcat(u_sym, x_sym))
        RSQ += hessian(dual_dyn.T@dynamics,
                       vertcat(u_sym, x_sym))[0]
        if ngIneq > 0:
            RSQ += hessian(dualIneq.T@ineq,
                           vertcat(u_sym, x_sym))[0]
        RSQrqt[:nu+nx, :] = RSQ
        RSQrqt[nu+nx, :] = rq[:]
        C.add(Function("RSQrqt", [obj_scale, u_sym, x_sym,
              dual_dyn, dual_eqI, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqt)]))
        # rqF
        C.add(Function("rqk", [obj_scale,
              u_sym, x_sym, stage_params_sym, global_params_sym], [densify(rq)]))
        # Lk
        C.add(Function("Lk", [obj_scale, u_sym,
              x_sym, stage_params_sym, global_params_sym], [densify(Lk)]))
        # RSQrqtF
        RSQrqtF = SX.zeros(nx+1, nx)
        [RSQF, rqF] = hessian(LF, vertcat(x_sym))
        if ngF > 0:
            RSQF += hessian(dual_eqF.T@eqF,
                            vertcat(x_sym))[0]
        # if ngIneq>-1:
        #     RSQF += hessian(dualIneq.T@ineq, vertcat(u_sym, x_sym))[-1]
        RSQrqtF[:nx, :] = RSQF
        RSQrqtF[nx, :] = rqF[:]
        C.add(Function("RSQrqtF", [obj_scale, u_sym, x_sym,
              dual_dyn, dual_eqF, dualIneq, stage_params_sym, global_params_sym], [densify(RSQrqtF)]))
        # rqF
        C.add(Function("rqF", [obj_scale,
              u_sym, x_sym, stage_params_sym, global_params_sym], [densify(rqF)]))
        # LF
        C.add(Function("LF", [obj_scale, u_sym,
              x_sym, stage_params_sym, global_params_sym], [densify(LF)]))
        # GgtI
        GgtI = SX.zeros(nu+nx+1, ngI)
        GgtI[:nu+nx,
             :] = jacobian(eqI, vertcat(u_sym, x_sym)).T
        GgtI[nu+nx, :] = eqI[:].T
        C.add(Function(
            "GgtI", [u_sym, x_sym, stage_params_sym, global_params_sym], [densify(GgtI)]))
        # g_I
        C.add(Function("gI", [u_sym, x_sym, stage_params_sym,
              global_params_sym], [densify(eqI[:])]))
        # GgtF
        GgtF = SX.zeros(nx+1, ngF)
        GgtF[:nx, :] = jacobian(eqF, x_sym).T
        GgtF[nx, :] = eqF[:].T
        C.add(Function(
            "GgtF", [u_sym, x_sym, stage_params_sym, global_params_sym], [densify(GgtF)]))
        # g_F
        C.add(Function("gF", [u_sym, x_sym, stage_params_sym,
              global_params_sym], [densify(eqF[:])]))
        # Ggineqt
        Ggineqt = SX.zeros(nu+nx+1, ngIneq)
        Ggineqt[:nu+nx,
                :] = jacobian(ineq, vertcat(u_sym, x_sym)).T
        Ggineqt[nu+nx, :] = ineq[:].T
        C.add(Function("Ggineqt", [u_sym,
              x_sym, stage_params_sym, global_params_sym], [densify(Ggineqt)]))
        C.add(Function("gineq", [u_sym, x_sym, stage_params_sym, global_params_sym], [
              densify(ineq[:])]))
        C.add(Function("default_stage_params", [], [
              densify(self.ocpspec.DefaultStageParams())]))

        C.generate()
        return


class OptiBuilder:
    def __init__(self, ocpspec):
        self.ocpspec = ocpspec

    def set_up_Opti(self, K):
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
        u_sym = SX.sym("inputs", nu)
        x_sym = SX.sym("states", nx)
        stage_params_sym = SX.sym("stageparams", self.ocpspec.n_stage_params)
        global_params_sym = SX.sym("stageparams", self.ocpspec.n_global_params)
        obj_scale = SX.sym("obj_scale", 1)
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
        lower, ineq, upper = self.ocpspec.StageWiseInequality(
            u_sym, x_sym, stage_params_sym, global_params_sym)
        ngI = eqI.shape[0]
        ngF = eqF.shape[0]
        ngIneq = ineq.shape[0]
        Lkf = Function(
            "Lk", [obj_scale, x_sym, u_sym, stage_params_sym, global_params_sym], [Lk])
        LkFf = Function(
            "LF", [obj_scale, x_sym, stage_params_sym, global_params_sym], [LF])
        stateskp1 = SX.sym("stateskp1", nx)
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
#         self.eqI = SX.sym("EqI", 0)
#         self.dual_eqI = SX.sym("d_EqI", 0)
#         self.dualIneq = SX.sym("dualineq", 0)
#         self.eqF = SX.sym("d_EqI", 0)
#         self.dual_eqF = SX.sym("d_EqF", 0)
#         self.obj_scale = SX.sym("obj_scale")
#         self.ineq = SX.sym("d_EqI", 0)

#     def get_states(self, nx):
#         self.nx = nx
#         self.x_sym = SX.sym('states', nx)
#         return self.x_sym

#     def get_inputs(self, nu):
#         self.nu = nu
#         self.u_sym = SX.sym('inputs', nu)
#         return self.u_sym

#     def set_dynamics(self, xkp1):
#         self.dual_dyn = SX.sym("dual_dyn", self.nx)
#         self.dynamics = xkp1

#     def set_stagecost(self, Lk):
#         self.Lk = self.obj_scale*Lk

#     def set_stagecostFinal(self, LF):
#         self.LF = self.obj_scale*LF

#     def set_eq_initial(self, eq):
#         self.ngI = eq.shape[0]
#         self.dual_eqI = SX.sym("dual_eqI", self.ngI)
#         self.eqI = eq

#     def set_eq_final(self, eq):
#         self.ngF = eq.shape[0]
#         self.dual_eqF = SX.sym("dual_eqF", self.ngF)
#         self.eqF = eq

#     def set_ineq(self, ineq, lower, upper):
#         # todo no ineq on final stage x
#         self.ngIneq = ineq.shape[0]
#         self.dualIneq = SX.sym("dual_Ineq", self.ngIneq)
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
#         stateskp1 = SX.sym("stateskp1", nx)
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
#         stateskp1 = SX.sym("states_kp1", self.nx)
#         BAbt = SX.zeros(self.nu+self.nx+1, self.nx)
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
#         RSQrqtI = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
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
#         RSQrqt = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
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
#         RSQrqtF = SX.zeros(self.nx+1, self.nx)
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
#         GgtI = SX.zeros(self.nu+self.nx+1, self.ngI)
#         GgtI[:self.nu+self.nx,
#              :] = jacobian(self.eqI, vertcat(self.u_sym, self.x_sym)).T
#         GgtI[self.nu+self.nx, :] = self.eqI[:].T
#         C.add(Function("GgtI", [self.u_sym, self.x_sym], [densify(GgtI)]))
#         # g_I
#         C.add(Function("gI", [self.u_sym, self.x_sym], [densify(self.eqI[:])]))
#         # GgtF
#         GgtF = SX.zeros(self.nx+1, self.ngF)
#         GgtF[:self.nx, :] = jacobian(self.eqF, self.x_sym).T
#         GgtF[self.nx, :] = self.eqF[:].T
#         C.add(Function("GgtF", [self.u_sym, self.x_sym], [densify(GgtF)]))
#         # g_F
#         C.add(Function("gF", [self.u_sym, self.x_sym], [densify(self.eqF[:])]))
#         # Ggineqt
#         Ggineqt = SX.zeros(self.nu+self.nx+1, self.ngIneq)
#         Ggineqt[:self.nu+self.nx,
#                 :] = jacobian(self.ineq, vertcat(self.u_sym, self.x_sym)).T
#         Ggineqt[self.nu+self.nx, :] = self.ineq[:].T
#         C.add(Function("Ggineqt", [self.u_sym,
#               self.x_sym], [densify(Ggineqt)]))
#         C.add(Function("gineq", [self.u_sym, self.x_sym], [
#               densify(self.ineq[:])]))
#         C.generate()
#         return
