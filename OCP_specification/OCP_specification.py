from casadi import *


class OptimalControlProblem:
    def __init__(self):
        self.ngI = 0
        self.ngF = 0
        self.obj_scale = SX.sym("obj_scale")

    def get_states(self, nx):
        self.nx = nx
        self.x_sym = SX.sym('states', nx)
        self.scales_x_sym = self.scales_x_sym*SX.sym('sclaes_states', nx)
        return self.x_sym

    def get_inputs(self, nu):
        self.nu = nu
        self.scales_u_sym = SX.sym('scales_states', nu)
        self.u_sym = self.scales_u_sym * SX.sym('states', nu)
        return self.u_sym

    def set_dynamics(self, xkp1):
        self.scales_dyn = SX.sym("scales_dyn", self.nx)
        self.dual_dyn = SX.sym("dual_dyn", self.nx)
        self.dynamics = self.scales_dyn * xkp1

    def set_stagecost(self, Lk):
        self.Lk = self.obj_scale*Lk

    def set_stagecostFinal(self, Lk):
        self.Lk = self.obj_scale*Lk

    def set_eq_initial(self, eq):
        self.ngI = eq.shape[0]
        self.scales_eqI = SX.sym("scales_eqI", self.ngI)
        self.dual_eqI = SX.sym("dual_eqI", self.ngI)
        self.eqI = self.scales_eqI * eq

    def set_eq_final(self, eq):
        self.ngF = eq.shape[0]
        self.scales_eqF = SX.sym("scales_eqI", self.ngF)
        self.dual_eqF = SX.sym("dual_eqF", self.ngF)
        self.eqF = self.scales_eqF*eq

    def generate_code(self, filename):
        C = CodeGenerator(filename)
        # BAbt
        stateskp1 = SX.sym("states_kp1", self.nx)
        scales_stateskp1 = SX.sym("scales_states_kp1", self.nx)
        BAbt = SX.zeros(self.nu+self.nx+1, self.nx)
        BAbt[:self.nu+self.nx,
             :] = jacobian(self.dynamics, vertcat(self.u_sym, self.x_sym)).T
        BAbt[self.nu+self.nx, :] = (stateskp1 - self.dynamics)[:]
        C.add(Function("BAbt", [stateskp1, scales_stateskp1, self.x_sym,
              self.scales_x_sym, self.u_sym, self.scales_u_sym, self.scales_dyn], [densify(BAbt)]))
        # b
        # RSQrqtI
        RSQrqtI = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQI, rqI] = hessian(self.Lk, vertcat(self.u_sym, self.x_sym))
        if self.ngI>0:
            RSQI += hessian(self.dual_eqI.T@self.eqI,
                            vertcat(self.u_sym, self.x_sym))
        RSQI += hessian(self.dual_dyn.T@self.dynamics,
                        vertcat(self.u_sym, self.x_sym))
        RSQrqtI[:self.nu+self.nx, :] = RSQI 
        RSQrqtI[self.nu+self.nx, :] = rqI[:] 
        C.add(Function("RSQrqtI", [self.obj_scale, self.x_sym, self.scales_x_sym, self.u_sym, self.scales_u_sym, self.dual_dyn, self.scales_dyn, self.dual_eqI, self.scales_eqI], [densify(RSQrqtI)]))
        # rqI
        # RSQrqt
        RSQrqt = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQ, rq] = hessian(self.Lk, vertcat(self.u_sym, self.x_sym))
        RSQ += hessian(self.dual_dyn.T@self.dynamics,
                        vertcat(self.u_sym, self.x_sym))
        RSQrqt[:self.nu+self.nx, :] = RSQ 
        RSQrqt[self.nu+self.nx, :] = rq[:] 
        C.add(Function("RSQrqt", [self.obj_scale, self.x_sym, self.scales_x_sym, self.u_sym, self.scales_u_sym, self.dual_dyn, self.scales_dyn, self.dual_eqI, self.scales_eqI], [densify(RSQrqt)]))
        # rq
        # RSQrqtF
        RSQrqtF = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQF, rqF] = hessian(self.Lk, vertcat(self.x_sym))
        if self.ngF>0:
            RSQF += hessian(self.dual_eqF.T@self.eqF,
                            vertcat(self.x_sym))
        RSQrqtF[self.nx:, :] = RSQF 
        RSQrqtF[self.nx, :] = rqF[:] 
        C.add(Function("RSQrqtF", [self.obj_scale, self.x_sym, self.scales_x_sym, self.u_sym, self.scales_u_sym, self.dual_dyn, self.scales_dyn, self.dual_eqF, self.scales_eqF], [densify(RSQrqtF)]))
        # rqF
        # GgtI
        GgtI = SX.zeros(self.nu+self.nx+1, self.ngI)
        GgtI[:self.nu+self.nx,:] = jacobian(self.eqI, vertcat(self.u_sym, self.x_sym)).T
        GgtI[self.nu+self.nx,:] = self.eqI[:]
        C.add(Function("GgtI", [self.x_sym, self.scales_x_sym, self.u_sym, self.scales_u_sym, self.scales_eqI], [densify(GgtI)]))
        # g_I
        # GgtF
        GgtF = SX.zeros(self.nx+1, self.ngI)
        GgtF[:self.nx,:] = jacobian(self.eqF, vertcat(self.x_sym)).T
        GgtF[self.nx,:] = self.eqF[:]
        C.add(Function("GgtI", [self.x_sym, self.scales_x_sym, self.u_sym, self.scales_u_sym, self.scales_eqI],[densify(GgtF)]))
        # g_F
        C.generate()
        pass
