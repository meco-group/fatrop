from casadi import *
import OCPSpecification
import RocketDynamics

OCP = OCPSpecification.OptimalControlProblem()
nu = 2
uk = OCP.get_inputs(nu)
xk = OCP.get_states(9)
robotdynamics = RocketDynamics.RockDyns()
OCP.set_dynamics(robotdynamics.dynamics(uk,xk)) #xkp1 = xk + uk
OCP.set_stagecost(xk.T@xk + sum1(xk) + uk.T@uk)
# OCP.set_stagecost( xk[0] )
OCP.set_stagecostFinal(0*xk.T@xk)
# OCP.set_eq_final(robotdynamics.xk_pos-DM.ones(2,1))
opti = OCP.set_up_Opti(10)
opti.solver("ipopt")
opti.solve()
opti.solve()
OCP.generate_code("f.c")
#gcc -fPIC -march=native -shared -O3 f.c -o f.so
