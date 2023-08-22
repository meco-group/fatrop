import fatropy
import numpy as np
rockit_generated_so = "foobar/casadi_codegen.so"
rockit_generated_json = "foobar/casadi_codegen.json"
ocp = fatropy.StageOCPApplicationFactory.from_rockit_interface(rockit_generated_so, rockit_generated_json)
ocp.optimize()
sol = ocp.last_solution()
# ocp.set_initial_x(sol.x)
# ocp.set_initial_u(sol.u)
ocp.set_initial(sol)
ocp.set_option("mu_init", 1e-3)
ocp.set_value("target", [1., 2])
# ocp.set_option("bound_push", 1e-6)
ocp.optimize()
print(ocp.last_solution().x)
# print(sol.evaluate(ocp.get_expression("control_u").at_control()))