from casadi import *
import OCPSpecification

OCP = OCPSpecification.OptimalControlProblem()
uk = OCP.get_inputs(2)
xk = OCP.get_states(3)
OCP.set_dynamics(vertcat(xk[:2]+uk, xk[2:])) #xkp1 = xk + uk
OCP.set_stagecost(xk.T@xk + sum1(xk) + 0.01*cos(uk.T@uk))
# OCP.set_stagecost( xk[0] )
OCP.set_stagecostFinal(xk.T@xk)
OCP.set_eq_final(xk[:2])
opti = OCP.set_up_Opti(10)
opti.solver("ipopt")
opti.solve()
OCP.generate_code("f.c")
#gcc -fPIC -march=native -shared -O3 f.c -o f.so
