from matplotlib.font_manager import json_load
from OCPSpecification import *
import TraGenSpec

import numpy as np
import json
FSGenData = json_load('FSGenData.json')
initial_frames = np.array(FSGenData["initial_frames"]).reshape(12,1824//12, order='c')
invariants = np.array(FSGenData["invariants_model"]).reshape(3, 152, order='c')
trajgenspec = TraGenSpec.TraGenSpec()
optibuilder = OptiBuilder(trajgenspec)
opti = optibuilder.set_up_Opti(152)
jsongen = JSONGenerator(trajgenspec)
R0 = initial_frames[:9, 0]
p0 = initial_frames[9:12, 0]
Rend = initial_frames[:9, -1]
pend = 1.5*initial_frames[9:12, -1]
stage_params = np.vstack((1.0/152.0*np.ones((1,152)),invariants))
global_params = np.hstack((R0,p0,Rend, pend))
initial_u = invariants[:,:-1]
initial_x = initial_frames
jsongen.generate_JSON("test.json", 152, stage_params, global_params, initial_x, initial_u)
opti.solver('ipopt',{}, {"print_level":5})
opti.set_initial(optibuilder.x_vars, initial_frames)
opti.set_initial(optibuilder.u_vars, invariants[:,:-1])
# todo: wrong should also include dt!
# print(np.vstack((1.0/152.0*np.ones((1,152)),invariants)))
opti.set_value(optibuilder.stage_params_in, np.vstack((1.0/152.0*np.ones((1,152)),invariants)))
opti.set_value(optibuilder.global_params_in, np.hstack((R0,p0,Rend, pend)))
opti.solve()
codegen = FatropOCPCodeGenerator(trajgenspec)
codegen.generate_code("f.c")