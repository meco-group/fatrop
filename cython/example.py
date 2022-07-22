import fatropy

functions = "../../fatrop-examples/rocket/Rocket_example.so"

json_spec = "../../fatrop-examples/rocket/Rocket_example.json"

myOCP = fatropy.PyOCP(functions,json_spec)

#print(myOCP.initial_x)
#print(myOCP.initial_u)
#myOCP.SetBounds()
#myOCP.SetInitial()

retval = myOCP.Optimize()
print(retval)

x_curr = myOCP.x_curr
x_next = myOCP.x_next

#print(x_curr)
#print(x_next)
print(myOCP.n_eqs)
print(myOCP.n_ineqs)
print(myOCP.total_time)