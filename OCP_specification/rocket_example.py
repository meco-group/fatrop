from casadi import *
import OCPSpecification

OCP = OCPSpecification.OptimalControlProblem()
uk = OCP.get_inputs(2)
xk = OCP.get_states(2)
OCP.set_dynamics(xk+uk) #xkp1 = xk + uk
OCP.set_stagecost(xk.T@xk + uk.T@uk)
OCP.set_stagecostFinal(xk.T@xk)
opti = OCP.set_up_Opti(10)
opti.solver("ipopt")
opti.solve()
OCP.generate_code("f.c")

