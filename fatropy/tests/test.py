import casadi as cs
import fatropy
x = cs.MX.sym("x")
func = cs.Function("test_func", [x], [x])
#fatropy.print_function_name(func)

