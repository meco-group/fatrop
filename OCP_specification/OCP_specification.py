from casadi import *

class OptimalControlProblem:
    def get_states(self, nx):
        self.nx = nx
        self.x_sym =  SX.sym('states', nx)
        return self.x_sym 
    def get_inputs(self, nu):
        self.nu = nu
        self.u_sym =  SX.sym('states', nu)
        return self.u_sym 
    def set_dynamics(self, xkp1):
        self.dynamics = Function('dyn_function', [self.u_sym,self.x_sym], [xkp1])
    def set_stagecost(self, Lk):
        self.Lk = Function('stagecost', [self.u_sym,self.x_sym], [Lk])
    def set_eq_initial(self, eq):
        self.eqI = Function('initial_c', [self.u_sym,self.x_sym], [eq])
    def set_eq_final(self, eq):
        self.eqF = Function('final_c', [self.x_sym], [eq])
    def generate_code(self, filename):
        C = CodeGenerator(filename)
        # BAbt
        # RSQrqI
        # RSQrqF
        # C.add(funci)
        C.generate()
        pass
