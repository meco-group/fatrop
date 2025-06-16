Fatrop CasADi Cookbook
======================

Fatrop is interfaced with CasADi's Opti framework, which provides a user-friendly way to define and solve optimization problems.
We'll walk through the process of setting up and solving a quadcopter trajectory optimization problem, focusing on key aspects and best practices.
This document is based on the example provided in the Fatrop repository, specifically `examples/ocp_quadrotor_example.py <https://github.com/meco-group/fatrop/blob/main/examples/ocp_quadrotor_example.py>`_.
It is not a complete copy of the example, but rather a guide to the key steps involved in using Fatrop with CasADi.

1. Problem Setup and Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we import the necessary libraries, define the physical parameters, and set up the dynamics of our quadcopter:

.. code-block:: python

   import casadi as ca
   import numpy as np

   # Physical parameters
   m = 1.0
   g = 9.81
   l = 0.25
   kF = 1e-5
   kM = 1e-6
   Ix, Iy, Iz = 0.01, 0.01, 0.02  # Diagonal inertia terms

   # rotor limits
   omega_min = 0.0
   omega_max = 600.0
   # orientation limits
   phi_min = -ca.pi / 4
   phi_max = ca.pi / 4
   theta_min = -ca.pi / 4
   theta_max = ca.pi / 4

   # State variables
   x, y, z = ca.MX.sym('x'), ca.MX.sym('y'), ca.MX.sym('z')
   vx, vy, vz = ca.MX.sym('vx'), ca.MX.sym('vy'), ca.MX.sym('vz')
   phi, theta, psi = ca.MX.sym('phi'), ca.MX.sym('theta'), ca.MX.sym('psi')
   p, q, r = ca.MX.sym('p'), ca.MX.sym('q'), ca.MX.sym('r')

   state = ca.vertcat(x, y, z, vx, vy, vz, phi, theta, psi, p, q, r)

   # Control inputs
   w1, w2, w3, w4 = ca.MX.sym('w1'), ca.MX.sym('w2'), ca.MX.sym('w3'), ca.MX.sym('w4')
   omega = ca.vertcat(w1, w2, w3, w4)

   # Define dynamics (state derivatives)
   state_dot = ...  # (full dynamics code here)

   # Create CasADi function for dynamics
   f_dyn = ca.Function('f', [state, omega], [state_dot]).expand()

   # Create CasADi function for integrator (explicit Euler)
   f_intg = ca.Function('f_dyn_discr', [state, omega], [state + T/K*f_dyn(state, omega)]).expand()

.. note::

   The integrator function will be re-used at every time step.  
   For functions which are used multiple times, it is generally beneficial to put them 
   in a CasADi function to avoid duplication of expressions.
   We use the expanded version of the dynamics function for efficiency.  
   Refer to the CasADi documentation for more details.

2. Constraints and Objective
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define constraints and the objective function:

.. code-block:: python

   # Control and state constraints
   omega_min, omega_max = 0.0, 1000.0
   phi_min, phi_max = -ca.pi/3, ca.pi/3
   theta_min, theta_max = -ca.pi/3, ca.pi/3

   # Initial and final states
   x0 = ca.DM([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
   xf = ca.DM([5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

   # Path constraints
   f_control = ca.Function('f_control', [state, omega], [ca.vertcat(omega, phi, theta)]).expand()
   lb_control = ca.vertcat(omega_min*ca.DM.ones(4), phi_min, theta_min)
   ub_control = ca.vertcat(omega_max*ca.DM.ones(4), phi_max, theta_max)

   # Objective function
   f_cost = ca.Function('f_cost', [omega], [ca.sumsqr(omega/omega_max)]).expand()

3. Problem Dimensions and Opti Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maintain a list of problem dimensions and set up the optimization problem.
For the CasADi's fatrop interface, it is important to define the optimization variables in a specific order:
x0, u0, x1, u1, ..., xK-1, uK-1.

.. code-block:: python

   K, T = 100, 5.0  # Number of control intervals, Total time
   nx = [12 for _ in range(K)] # number of state variables at each time step
   nu = [4 for _ in range(K-1)] + [0] # number of control inputs at each time step
   ng = []  # number of path constraints, rill we populated when setting up constraints

   opti = ca.Opti()
   x = []
   u = []
   for k in range(K):
       x.append(opti.variable(nx[k]))
       u.append(opti.variable(nu[k]))


4. Constraints, Objective, and Initial Guess
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set up constraints, objective, and provide an initial guess.
The constraints should be defined in the following order:
discrete_dynamics_0, path_constraints_0, ..., discrete_dynamics_K-2, path_constraints_K-2, path_constraints_K-1

.. code-block:: python

   ng = []
   for k in range(K):
       if k < K-1:
           opti.subject_to(x[k+1] == discrete_dynamics(u[k], x[k], k))
       path_constr = path_constraints(u[k], x[k], k)
       ng.append(0)
       for constr in path_constr:
           ng[-1] += constr[1].nnz()
           opti.subject_to((constr[0] <= constr[1]) <= constr[2])

   # Objective
   J = 0
   for k in range(K):
       J += cost(u[k], x[k], k)
   opti.minimize(J)

   # Initial guess
   for k in range(K):
       u_init, x_init = initial_guess(k)
       opti.set_initial(u[k], u_init)
       opti.set_initial(x[k], x_init)

The `ng` list keeps track of the number of path constraints at each time step, which is important for the manual structure detection when using Fatrop.

5. Solving with Ipopt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve the problem using both Ipopt (reference solve):

.. code-block:: python

   # Solve with Ipopt
   opti.solver('ipopt', {})
   sol_ipopt = opti.solve()

5. Solving with Fatrop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Solve the problem using both Ipopt and Fatrop:

.. code-block:: python

   # Solve with Fatrop
   opti.solver('fatrop', {
       'structure_detection': 'manual',
       'nx': nx, 'nu': nu, 'ng': ng, 'N': K-1,
       "expand": False})
   sol_fatrop = opti.solve()

6. Results and Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Retrieve and visualize the results:

.. code-block:: python

   X_sol = sol_fatrop.value(ca.horzcat(*x))

   # Visualization code
   # ... (3D trajectory plotting)

For the complete visualization code, please refer to the full example in `examples/ocp_quadrotor_example.py <https://github.com/meco-group/fatrop/blob/main/examples/ocp_quadrotor_example.py>`_.

Advanced Usage and Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expanding Functions
-------------------

To potentially speed up function evaluation, we've used expanded functions throughout this example:

.. code-block:: python

   f_dyn = ca.Function('f', [state, omega], [state_dot]).expand()
   f_cost = ca.Function('f_cost', [omega], [ca.sum1((omega/omega_max)**2)]).expand()
   f_control = ca.Function('f_control', [state, omega], [ca.vertcat(omega, phi, theta)]).expand()

Performance might also improve by expanding the full functions, used by fatrop internally.
This can be done by setting the `expand` option to `True` when creating the function.

.. code-block:: python

   # Solve with Fatrop
   opti.solver('fatrop', {
       'structure_detection': 'manual',
       'nx': nx, 'nu': nu, 'ng': ng, 'N': K-1,
       "expand": True})
   sol_fatrop = opti.solve()

This can speed up the function evaluation, sometimes at the cost of having larger expressions with duplicated code.
These large expressions can result in long compilation times when using code generation / JIT compilation.

Just-in Time (JIT) Compilation of Functions 
-------------------------------------------

CasADi supports Just-in-Time (JIT) compilation, which can significantly speed up the evaluation of functions.

.. code-block:: python

   opti.solver('fatrop', {'structure_detection':'manual', 'nx': nx, 'nu':nu, 'ng':ng, 'N':K-1, "expand": False, "jit":True, "jit_options": {"flags": "-O3", "verbose": True}})
   res = opti.solve()

Code Generation
---------------

Code generation is a powerful technique to improve the performance of your optimization problem.
CasADi provides tools to generate C code for any CasADi function.
This means that we can put the full solver in a CasADi function and generate C code for it.

1. Generate C code for the optimization problem:

   .. code-block:: python

      # Generate C code
      opti.to_function('f_quadrotor', [ca.horzcat(*x), ca.horzcat(*u[:-1])], [ca.horzcat(*x)], ['x', 'u'], ['X']).generate('quadrotor_ocp.c', {"with_header": True})

   This generates a C file and header named 'quadrotor_ocp.(c/h)' that contains the optimized code/declarations for your problem.

The generated code can be used directly from C/C++ applications, and it's completely CasADi-free.
This means you can integrate the optimized function into your C/C++ projects without needing CasADi as a dependency.
The generated code can be compiled with:

   .. code-block:: bash

      gcc -fPIC -shared quadcopter.c -g -O3 -march=native -lfatrop -lblasfeo -I`fatrop path` -I`blasfeo path`/include/blasfeo/include 

This shared library can be imported into casadi using CasADi's `external` function interface.

For more information refer to the directory `examples/casadi_codegen/ <https://github.com/meco-group/fatrop/tree/main/examples/casadi_codegen>`_ in the Fatrop repository.


Further References
--------------------------------

For more information on using Fatrop with CasADi, refer to the following resources:
 - `Fatrop CasADi video tutorial on YouTube <https://www.youtube.com/watch?v=example>`_
 - `Fatrop CasADi demo Github repo <https://github.com/jgillis/fatrop_demo>`_
 - `CasADi website <https://web.casadi.org/>`_