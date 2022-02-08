from casadi import *


# This class can be used to specify OCP's with intial and final constraints and generate necessary code for it.
class OptimalControlProblem:
    def __init__(self):
        self.ngI = 0
        self.ngF = 0
        self.eqI = SX.sym("EqI", 0)
        self.dual_eqI = SX.sym("d_EqI", 0)
        self.eqF = SX.sym("d_EqI", 0)
        self.dual_eqF = SX.sym("d_EqF", 0)
        self.obj_scale = SX.sym("obj_scale")

    def get_states(self, nx):
        self.nx = nx
        self.x_sym = SX.sym('states', nx)
        return self.x_sym

    def get_inputs(self, nu):
        self.nu = nu
        self.u_sym = SX.sym('inputs', nu)
        return self.u_sym

    def set_dynamics(self, xkp1):
        self.dual_dyn = SX.sym("dual_dyn", self.nx)
        self.dynamics = xkp1

    def set_stagecost(self, Lk):
        self.Lk = self.obj_scale*Lk

    def set_stagecostFinal(self, LF):
        self.LF = self.obj_scale*LF

    def set_eq_initial(self, eq):
        self.ngI = eq.shape[0]
        self.dual_eqI = SX.sym("dual_eqI", self.ngI)
        self.eqI = eq

    def set_eq_final(self, eq):
        self.ngF = eq.shape[0]
        print(self.ngF)
        self.dual_eqF = SX.sym("dual_eqF", self.ngF)
        self.eqF = eq

    def set_up_Opti(self, K):
        # all scales are set to 1.0
        nu = self.nu
        nx = self.nx
        self.opti = Opti()
        self.N_vars = K*nx + (K-1)*nu
        self.opti_vars = self.opti.variable(self.N_vars)
        self.x_vars = MX.zeros(nx, K)
        self.u_vars = MX.zeros(nu, K)
        for k in range(K-1):
            offs = k*(nu+nx)
            self.u_vars[:, k] = self.opti_vars[offs:offs+nu]
            self.x_vars[:, k] = self.opti_vars[offs+nu:offs+nu+nx]
        offs = (K-1)*(nu+nx)
        self.x_vars[:, K-1] = self.opti_vars[offs:offs+nx]
        Lkf = Function(
            "Lk", [self.obj_scale, self.x_sym, self.u_sym], [self.Lk])
        LkFf = Function("LF", [self.obj_scale, self.x_sym], [self.LF])
        stateskp1 = SX.sym("stateskp1", nx)
        Dynamcisf = Function(
            "F", [stateskp1, self.x_sym, self.u_sym], [-stateskp1 + self.dynamics])
        if self.ngI > 0:
            EqIf = Function("eqI", [self.x_sym, self.u_sym], [self.eqI])
            self.opti.subject_to(
                EqIf(self.x_vars[:, 0], self.u_vars[:, 0]) == 0.0)
        if self.ngF > 0:
            EqFf = Function("eqF", [self.x_sym], [self.eqF])
            self.opti.subject_to(EqFf(self.x_vars[:, K-1]) == 0.0)
        J = 0
        for k in range(K-1):
            J += Lkf(1.0, self.x_vars[:, k], self.u_vars[:, k])
            self.opti.subject_to(
                Dynamcisf(self.x_vars[:, k+1], self.x_vars[:, k], self.u_vars[:, k]) == 0.0)
        J += LkFf(1.0, self.x_vars[:, K-1])
        self.opti.minimize(J)
        return self.opti

    def generate_code(self, filename):
        C = CodeGenerator(filename)
        # BAbt
        stateskp1 = SX.sym("states_kp1", self.nx)
        BAbt = SX.zeros(self.nu+self.nx+1, self.nx)
        BAbt[:self.nu+self.nx,
             :] = jacobian(self.dynamics, vertcat(self.u_sym, self.x_sym)).T
        b = (-stateskp1 +  self.dynamics)[:]
        BAbt[self.nu+self.nx, :] = b
        C.add(
            Function("BAbt", [stateskp1, self.u_sym, self.x_sym], [densify(BAbt)]))
        # b
        C.add(Function("bk", [stateskp1, self.u_sym,
                              self.x_sym], [densify(b)]))
        # RSQrqtI
        RSQrqtI = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQI, rqI] = hessian(self.Lk, vertcat(self.u_sym, self.x_sym))
        if self.ngI > 0:
            RSQI += hessian(self.dual_eqI.T@self.eqI,
                            vertcat(self.u_sym, self.x_sym))[0]
        RSQI += hessian(self.dual_dyn.T@self.dynamics,
                        vertcat(self.u_sym, self.x_sym))[0]
        RSQrqtI[:self.nu+self.nx, :] = RSQI
        RSQrqtI[self.nu+self.nx, :] = rqI[:]
        C.add(Function("RSQrqtI", [self.obj_scale, self.u_sym,
              self.x_sym, self.dual_dyn, self.dual_eqI], [densify(RSQrqtI)]))
        rqI
        C.add(Function("rqI", [self.obj_scale,
              self.u_sym, self.x_sym], [densify(rqI)]))
        # RSQrqt
        RSQrqt = SX.zeros(self.nu+self.nx+1, self.nu + self.nx)
        [RSQ, rq] = hessian(self.Lk, vertcat(self.u_sym, self.x_sym))
        RSQ += hessian(self.dual_dyn.T@self.dynamics,
                       vertcat(self.u_sym, self.x_sym))[0]
        RSQrqt[:self.nu+self.nx, :] = RSQ
        RSQrqt[self.nu+self.nx, :] = rq[:]
        C.add(Function("RSQrqt", [self.obj_scale, self.u_sym, self.x_sym,self.dual_dyn, self.dual_eqI], [densify(RSQrqt)]))
        # rqF
        C.add(Function("rqk", [self.obj_scale, self.u_sym, self.x_sym], [densify(rq)]))
        # Lk
        C.add(Function("Lk", [self.obj_scale, self.u_sym, self.x_sym], [densify(self.Lk)]))
        # RSQrqtF
        RSQrqtF = SX.zeros(self.nx+1, self.nx)
        [RSQF, rqF] = hessian(self.LF, vertcat(self.x_sym))
        if self.ngF > 0:
            RSQF += hessian(self.dual_eqF.T@self.eqF,
                            vertcat(self.x_sym))[0]
        RSQrqtF[:self.nx, :] = RSQF
        RSQrqtF[self.nx, :] = rqF[:]
        C.add(Function("RSQrqtF", [self.obj_scale, self.u_sym, self.x_sym,
              self.dual_dyn, self.dual_eqF], [densify(RSQrqtF)]))
        # rqF
        C.add(Function("rqF", [self.obj_scale, self.u_sym, self.x_sym], [densify(rqF)]))
        # LF
        C.add(Function("LF", [self.obj_scale, self.u_sym, self.x_sym], [densify(self.LF)]))
        # GgtI
        GgtI = SX.zeros(self.nu+self.nx+1, self.ngI)
        GgtI[:self.nu+self.nx,
             :] = jacobian(self.eqI, vertcat(self.u_sym, self.x_sym)).T
        GgtI[self.nu+self.nx, :] = self.eqI[:].T
        C.add(Function("GgtI", [self.u_sym, self.x_sym], [densify(GgtI)]))
        # g_I
        C.add(Function("gI", [self.u_sym, self.x_sym], [densify(self.eqI[:])]))
        # GgtF
        GgtF = SX.zeros(self.nx+1, self.ngF)
        GgtF[:self.nx, :] = jacobian(self.eqF, self.x_sym).T
        GgtF[self.nx, :] = self.eqF[:].T
        C.add(Function("GgtF", [self.u_sym, self.x_sym], [densify(GgtF)]))
        # g_F
        C.add(Function("gF", [self.u_sym, self.x_sym], [densify(self.eqF[:])]))
        C.generate()
        return
