from fatropy.spectool import *
from casadi import *
def get_func():
  ocp = Ocp()
  x = ocp.state(2)
  u = ocp.control()
  p = ocp.parameter()

  e = 1. - x[0] * x[0]
  dt = .5
  x_next = x + vertcat(x[1], e * x[1] - x[0] + u) * dt

  # 
  #   =----  initial stage ----=
  # 
  initial_ustage = ocp.new_ustage() # states and controls are derived automatically
  # constraint 
  initial_ustage.subject_to(x[0] == 1.)
  # constraint 
  initial_ustage.subject_to(x[1] == 0.)
  # dynamics 
  initial_ustage.set_next(x, x_next)
  # objective  
  initial_ustage.add_objective(u * u + sumsqr(x))

  # 
  #   =----  middle stage ----=
  # 
  middle_ustage = ocp.new_ustage(19) # 19 is the number of nodes, states, controls and parameters are derived automatically
  # constraints 
  middle_ustage.subject_to(-0.25 < x[1])
  middle_ustage.subject_to((-1.0 < u) < 1)
  # dynamics  
  middle_ustage.set_next(x, x_next)
  # objective 
  middle_ustage.add_objective(u * u + sumsqr(x))

  # 
  #   =----  terminal stage ----=
  # 
  terminal_ustage = ocp.new_ustage() # the last stage also has x2 as state
  # constraints  
  terminal_ustage.subject_to(-0.25 < x[1])
  terminal_ustage.subject_to(x[1] == p)
  # objective 
  terminal_ustage.add_objective(x[1]*x[1]+p)

  # ocp.set_initial(u, p)

  return ocp.to_function([p], [ocp.at_t0(u), ocp.sample(x)[0],p], {"jit":True})

ocp_func = get_func()
ret = ocp_func(2.5)
print(ret)