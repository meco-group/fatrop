import fatropy
import numpy as np
ocp = fatropy.StageOCPApplicationFactory.from_rockit_interface("/home/lander/projects/fatrop_rockit_demo/foobar/casadi_codegen.so", "/home/lander/projects/fatrop_rockit_demo/foobar/casadi_codegen.json")
ocp.optimize()
sol = ocp.last_stageocp_solution()
# ocp.set_initial_x(sol.x)
# ocp.set_initial_u(sol.u)
ocp.set_initial(sol)
ocp.set_option("mu_init", 1e-3)
ocp.set_value("target", [1., 2])
# ocp.set_option("bound_push", 1e-6)
ocp.optimize()
print(ocp.last_stageocp_solution().x)
# print(sol.evaluate(ocp.get_expression("control_u").at_control()))