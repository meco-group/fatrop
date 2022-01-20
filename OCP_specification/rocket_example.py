from casadi import *
import OCPSpecification

OCP = OCPSpecification.OptimalControlProblem()
uk = OCP.get_inputs(3)
xk = OCP.get_states(12)
OCP.set_dynamics(vertcat(xk[:3]+uk, xk[3:])) #xkp1 = xk + uk
OCP.set_stagecost(xk.T@xk + uk.T@uk)
OCP.set_stagecostFinal(xk.T@xk)
opti = OCP.set_up_Opti(10)
opti.solver("ipopt")
opti.solve()
OCP.generate_code("f.c")
#gcc -fPIC -march=native -shared -O3 f.c -o f.so
