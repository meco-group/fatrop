from casadi import *
import FatropOCPSpecification
import RocketDynamics
import RocketSpec

rocketspec = RocketSpec.RocketSpec()
codegen = FatropOCPSpecification.FatropOCPCodeGenerator(rocketspec)
codegen.generate_code("f.c")
optibuilder = FatropOCPSpecification.OptiBuilder(rocketspec)
jsongen = FatropOCPSpecification.JSONGenerator(rocketspec)
K = 100
jsongen.generate_JSON("rocket_example.json", K, np.zeros((0,K)), np.zeros((0,K)), np.ones((rocketspec.nx, K)), np.ones((rocketspec.nu,K-1)))
opti = optibuilder.set_up_Opti(K)
opti.set_initial(optibuilder.u_vars, DM.ones(rocketspec.nu, K-1))
opti.set_initial(optibuilder.x_vars, DM.ones(rocketspec.nx, K))
opti.solver("ipopt")
opti.solve()

# OCP = OCPSpecification.OptimalControlProblem()
# nu = 2
# uk = OCP.get_inputs(nu)
# xk = OCP.get_states(9)
# robotdynamics = RocketDynamics.RockDyns()
# OCP.set_dynamics(robotdynamics.dynamics(uk,xk)) #xkp1 = xk + uk
# OCP.set_stagecost(xk.T@xk + sum1(xk) + uk.T@uk )
# # OCP.set_stagecost( xk[0] )
# OCP.set_stagecostFinal(xk.T@xk)
# OCP.set_eq_final(vertcat(robotdynamics.xk_pos-10*DM.ones(2,1),robotdynamics.xk_vel))
# OCP.set_eq_initial(vertcat(robotdynamics.xk_pos,robotdynamics.xk_vel))
# OCP.set_ineq(uk[0], [0], [inf])
# # OCP.set_ineq(uk[0], [-inf], [0])
# opti = OCP.set_up_Opti(100)
# opti.solver("ipopt")
# # opti.solver("ipopt",{}, {"print_level":5, "max_iter":500, "max_soc":5, "constr_mult_init_max":0.0, "bound_relax_factor":0.0, "kappa_d":0.0})
# # opti.solver("ipopt",{}, {"print_level":5, "max_iter":500, "max_soc":0,  "bound_relax_factor":0.0, "kappa_d":0.0})
# opti.solver("ipopt",{}, {"print_level":5, "max_iter":500, "max_soc":0,  "bound_relax_factor":0.0, "kappa_d":1e-4, "print_timing_statistics":"yes"})
# # opti.solver("ipopt",{}, {"print_level":7, "max_iter":3, "max_soc":0,  "bound_relax_factor":0.0, "kappa_d":1e-4})
# opti.set_initial(OCP.opti_vars, DM.ones(1,OCP.N_vars))
# opti.solve()
# OCP.generate_code("f.c")
# gcc -fPIC -march=native -shared -O3 f.c -o f.so
