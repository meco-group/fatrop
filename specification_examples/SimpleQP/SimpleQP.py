from QPSpecification import *
import numpy as np
import numpy.matlib
from FatropOCPSpecification import *
K = 20 
qpspecification = QPSpecification()
optibuilder = OptiBuilder(qpspecification)
opti = optibuilder.set_up_Opti(K)
stage_params = np.zeros((0,20))
global_params = np.array([])
inits_x = np.zeros((qpspecification.nx, K))
inits_u = np.zeros((qpspecification.nu, K-1))
opti.set_value(optibuilder.stage_params_in, stage_params)
opti.set_value(optibuilder.global_params_in, global_params)
opti.set_initial(optibuilder.x_vars, inits_x)
opti.set_initial(optibuilder.u_vars, inits_u)
opti.solver("ipopt", {"expand":True}, {'print_level':5, 'tol':1e-5, 'max_iter':500, "max_soc":0, "min_refinement_steps":0, 'bound_relax_factor':0.0, 'nlp_scaling_method':'none', 'kappa_sigma':1e10})
opti.solve()
# opti.solve()
lowerF = qpspecification.lowerF
upperF = qpspecification.upperF
jsongen = JSONGenerator(qpspecification)
jsongen.generate_JSON('test.json', K, stage_params, global_params, inits_x, inits_u, np.zeros((0, K)), np.zeros((0,K)), lowerF, upperF)
codegen = FatropOCPCodeGenerator(qpspecification)
codegen.generate_code("f.c")