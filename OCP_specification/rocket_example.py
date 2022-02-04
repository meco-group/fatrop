from casadi import *
import OCPSpecification
import RocketDynamics

OCP = OCPSpecification.OptimalControlProblem()
nu = 2
uk = OCP.get_inputs(nu)
xk = OCP.get_states(9)
robotdynamics = RocketDynamics.RockDyns()
OCP.set_dynamics(robotdynamics.dynamics(uk,xk)) #xkp1 = xk + uk
OCP.set_stagecost(xk.T@xk + sum1(xk) + uk.T@uk )
# OCP.set_stagecost( xk[0] )
OCP.set_stagecostFinal(xk.T@xk)
<<<<<<< Updated upstream
# OCP.set_eq_final(vertcat(robotdynamics.xk_pos-11*DM.ones(2,1),robotdynamics.xk_vel))
# OCP.set_eq_initial(vertcat(robotdynamics.xk_pos,robotdynamics.xk_vel))
opti = OCP.set_up_Opti(100)
opti.solver("ipopt")
# opti.solver("ipopt",{}, {"print_level":7, "max_iter":3, "max_soc":1})
=======
OCP.set_eq_final(vertcat(robotdynamics.xk_pos-10*DM.ones(2,1),robotdynamics.xk_vel))
# OCP.set_eq_initial(vertcat(robotdynamics.xk_pos,robotdynamics.xk_vel))
opti = OCP.set_up_Opti(100)
opti.solver("ipopt",{}, {"print_level":5, "max_iter":100, "max_soc":0})
>>>>>>> Stashed changes
opti.set_initial(OCP.opti_vars, DM.ones(1,OCP.N_vars))
opti.solve()
OCP.generate_code("f.c")
#gcc -fPIC -march=native -shared -O3 f.c -o f.so
