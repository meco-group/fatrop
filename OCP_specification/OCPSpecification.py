from casadi import *


# This class can be used to specify OCP's with intial and final constraints and generate necessary code for it.
class OptimalControlProblem:
    def __init__(self):
        self.ngI = 0
        self.ngF = 0
        self.eqI = SX.sym("EqI", 0)
        self.dual_eqI = SX.sym("d_EqI", 0)
        self.scales_eqI = SX.sym("d_EqI", 0)
        self.eqF = SX.sym("d_EqI", 0)
        self.dual_eqF = SX.sym("d_EqF", 0)
        self.scales_eqF = SX.sym("d_EqF", 0)
        self.obj_scale = SX.sym("obj_scale")

    def get_states(self, nx):
        self.nx = nx
        self.scales_x_sym = SX.sym('scales_states', nx)
        self.x_sym_scaled = SX.sym('states', nx)
        self.x_sym = self.scales_x_sym*self.x_sym_scaled
        return self.x_sym

    def get_inputs(self, nu):
        self.nu = nu
        self.scales_u_sym = SX.sym('scales_inputs', nu)
        self.u_sym_scaled =  SX.sym('inputs', nu)
        self.u_sym = self.scales_u_sym * self.u_sym_scaled
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

    def set_up_Opti(self, K):
        # all scales are set to 1.0
        nu = self.nu
        nx = self.nx
        self.opti = Opti()
        N_vars = K*nx + (K-1)*nu
        self.opti_vars = self.opti.variable(N_vars)
        self.x_vars = MX.zeros(nx, K)
        self.u_vars = MX.zeros(nu, K)
        for k in range(K-1):
            offs = k*(nu+nx)
            self.u_vars[:, k] = self.opti_vars[offs:offs+nu]
            self.x_vars[:, k] = self.opti_vars[offs+nu:offs+nu+nx]
        offs = (K-1)*(nu+nx)
        self.x_vars[:, K-1] = self.opti_vars[offs:offs+nx]
        Lkf = Function("Lk", [self.obj_scale, self.x_sym_scaled,
                       self.scales_x_sym, self.u_sym_scaled, self.scales_u_sym], [self.Lk])
        LkFf = Function("Lk", [self.obj_scale, self.x_sym_scaled,
                        self.scales_x_sym], [self.Lk])
        stateskp1 = SX.sym("stateskp1", nx)
        Dynamcisf = Function("F", [stateskp1, self.scales_x_sym, self.x_sym_scaled, self.u_sym_scaled,
                             self.scales_u_sym, self.scales_dyn], [-stateskp1 + self.dynamics])
        if self.ngI > 0:
            EqIf = Function("eqI", [self.x_sym_scaled, self.scales_x_sym,
                            self.u_sym_scaled, self.scales_u_sym, self.scales_eqI], [self.eqI])
            self.opti.subject_to(EqIf(self.x_vars[:, 0], DM.ones(
                nx), self.u_vars[:, 0], DM.ones(nu), DM.ones(self.ngI)) == 0.0)
        if self.ngF > 0:
            EqFf = Function(
                "eqI", [self.x_sym_scaled, self.scales_x_sym, self.scales_eqF], [self.eqF])
            self.opti.subject_to(
                EqFf(self.x_vars[:, K-1], DM.ones(nx), DM.ones(self.ngI)) == 0.0)
        J = 0
        for k in range(K-1):
            J += Lkf(1.0, self.x_vars[:, k],
                          DM.ones(nx), self.u_vars[:, k], DM.ones(nu))
            self.opti.subject_to(Dynamcisf(self.x_vars[:, k+1], DM.ones(
                nx), self.x_vars[:, k], self.u_vars[:, k], DM.ones(nu), DM.ones(nx))==0.0)
        J += LkFf(1.0, self.x_vars[:, K-1], DM.ones(nx))
        self.opti.minimize(J)
        return self.opti

    def generate_code(self, filename):
        C = CodeGenerator(filename)
        # BAbt
        stateskp1 = SX.sym("states_kp1", self.nx)
        scales_stateskp1 = SX.sym("scales_states_kp1", self.nx)
        BAbt = SX.zeros(self.nu+self.nx+1, self.nx)
        BAbt[:self.nu+self.nx,
             :] = jacobian(self.dynamics, vertcat(self.u_sym_scaled, self.x_sym_scaled)).T
        BAbt[self.nu+self.nx, :] = (stateskp1 - self.dynamics)[:]
        C.add(Function("BAbt", [stateskp1, scales_stateskp1, self.x_sym_scaled,
              self.scales_x_sym, self.u_sym_scaled, self.scales_u_sym, self.scales_dyn], [densify(BAbt)]))
        # b
        # RSQrqtI
        RSQrqtI = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQI, rqI] = hessian(self.Lk, vertcat(self.u_sym_scaled, self.x_sym_scaled))
        if self.ngI > 0:
            RSQI += hessian(self.dual_eqI.T@self.eqI,
                            vertcat(self.u_sym_scaled, self.x_sym_scaled))
        RSQI += hessian(self.dual_dyn.T@self.dynamics,
                        vertcat(self.u_sym_scaled, self.x_sym_scaled))[0]
        RSQrqtI[:self.nu+self.nx, :] = RSQI
        RSQrqtI[self.nu+self.nx, :] = rqI[:]
        C.add(Function("RSQrqtI", [self.obj_scale, self.x_sym_scaled, self.scales_x_sym, self.u_sym_scaled, self.scales_u_sym,
              self.dual_dyn, self.scales_dyn, self.dual_eqI, self.scales_eqI], [densify(RSQrqtI)]))
        # rqI
        # RSQrqt
        RSQrqt = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQ, rq] = hessian(self.Lk, vertcat(self.u_sym_scaled, self.x_sym_scaled))
        RSQ += hessian(self.dual_dyn.T@self.dynamics,
                       vertcat(self.u_sym_scaled, self.x_sym_scaled))[0]
        RSQrqt[:self.nu+self.nx, :] = RSQ
        RSQrqt[self.nu+self.nx, :] = rq[:]
        C.add(Function("RSQrqt", [self.obj_scale, self.x_sym_scaled, self.scales_x_sym, self.u_sym_scaled, self.scales_u_sym,
              self.dual_dyn, self.scales_dyn, self.dual_eqI, self.scales_eqI], [densify(RSQrqt)]))
        # rq
        # RSQrqtF
        RSQrqtF = SX.zeros(self.nx+1, self.nx)
        [RSQF, rqF] = hessian(self.Lk, vertcat(self.x_sym_scaled))
        if self.ngF > 0:
            RSQF += hessian(self.dual_eqF.T@self.eqF,
                            vertcat(self.x_sym_scaled))
        RSQrqtF[:self.nx, :] = RSQF
        RSQrqtF[self.nx, :] = rqF[:]
        C.add(Function("RSQrqtF", [self.obj_scale, self.x_sym_scaled, self.scales_x_sym, self.u_sym_scaled, self.scales_u_sym,
              self.dual_dyn, self.scales_dyn, self.dual_eqF, self.scales_eqF], [densify(RSQrqtF)]))
        # rqF
        # GgtI
        GgtI = SX.zeros(self.nu+self.nx+1, self.ngI)
        GgtI[:self.nu+self.nx,
             :] = jacobian(self.eqI, vertcat(self.u_sym_scaled, self.x_sym_scaled)).T
        GgtI[self.nu+self.nx, :] = self.eqI[:].T
        C.add(Function("GgtI", [self.x_sym_scaled, self.scales_x_sym, self.u_sym_scaled,
              self.scales_u_sym, self.scales_eqI], [densify(GgtI)]))
        # g_I
        # GgtF
        GgtF = SX.zeros(self.nx+1, self.ngI)
        GgtF[:self.nx, :] = jacobian(self.eqF, vertcat(self.x_sym_scaled)).T
        GgtF[self.nx, :] = self.eqF[:].T
        C.add(Function("GgtF", [self.x_sym_scaled, self.scales_x_sym, self.u_sym_scaled,
              self.scales_u_sym, self.scales_eqI], [densify(GgtF)]))
        # g_F
        C.generate()
        return
