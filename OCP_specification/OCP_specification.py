from casadi import *

class OptimalControlProblem:
    def get_states(self, nx):
        self.nx = nx
        self.x_sym =  SX.sym('states', nx)
        self.scales_x_sym = self.scales_x_sym*SX.sym('sclaes_states', nx)
        return self.x_sym 
    def get_inputs(self, nu):
        self.nu = nu
        self.scales_u_sym =  SX.sym('scales_states', nu)
        self.u_sym =  self.scales_u_sym * SX.sym('states', nu)
        return self.u_sym 
    def set_dynamics(self, xkp1):
        self.scales_dyn = SX.sym("scales_dyn", xkp1.shape[0])
        self.dynamics = self.scales_dyn * xkp1
    def set_stagecost(self, Lk):
        self.obj_scale = SX.sym("obj_scale")
        self.Lk = self.obj_scale*Lk
    def set_eq_initial(self, eq):
        self.ngI = eq.shape[0]
        self.scales_eqI = SX.sym("scales_eqI", self.ngI) 
        self.eqI = self.scales_eqI *eq
    def set_eq_final(self, eq):
        self.ngF = eq.shape[0]
        self.scales_eqF = SX.sym("scales_eqI", self.ngF) 
        self.eqF = self.scales_eqF*eq
    def generate_code(self, filename):
        C = CodeGenerator(filename)
        # BAbt
        stateskp1 = SX.sym("states_kp1", self.nx)
        scales_stateskp1 = SX.sym("scales_states_kp1", self.nx)
        BAbt = SX.zeros(self.nu+self.nx+1, self.nu+ self.nx)
        BAbt[:self.nu+self.nx, :] = jacobian(self.dynamics, vertcat(self.u_sym, self.x_sym)).T
        BAbt[self.nu+self.nx, :] = (stateskp1 - self.dynamics)[:]
        C.add(Function("BAbt", [stateskp1, scales_stateskp1, self.x_sym, self.scales_x_sym, self.u_sym, self.scales_u_sym, self.scales_dyn], [BAbt]))
        # RSQrqtI
        # rqI
        # RSQrqt
        # rq
        # RSQrqtF
        # rqF
        # GgtI 
        # eq_I 
        # GgtF 
        # eq_F 
        # C.add(funci)
        C.generate()
        pass
