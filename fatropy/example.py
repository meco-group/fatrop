import fatropy
import numpy as np
ocp = fatropy.StageOCPApplicationFactory.from_rockit_interface("/home/lander/projects/fatrop_rockit_demo/foobar/casadi_codegen.so", "/home/lander/projects/fatrop_rockit_demo/foobar/casadi_codegen.json")
print(ocp.optimize())
sol = ocp.last_solution()
ocp.set_initial(sol)
ocp.set_option("mu_init", 1e-3)
# ocp.set_value("target", [1., 2])
ocp.set_option("bound_push", 1e-6)
ocp.optimize()
x = np.array(ocp.last_solution().x)
print(x)
# print(sol.evaluate(ocp.get_expression("control_u").at_control()))