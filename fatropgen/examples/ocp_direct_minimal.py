"""Minimal OcpProxy example: a 1-D double integrator.

The smallest useful demonstration of specifying a Fatrop problem with the
high-level :class:`OcpProxy` and solving it via the direct-``OcpAbstract``
codegen.  Read this first to learn the proxy API.

Problem: move a point on a line from rest at ``p=0`` to rest at ``p=1`` while
minimising the sum of squared accelerations.

    state    x = [p, v]          control  u = [a]
    dynamics p+ = p + dt*v ,      v+ = v + dt*a
    stage cost   a^2
    start    x = [0, 0]   goal  x = [1, 0]   bound  -5 <= a <= 5

A proxy is just a list of *node templates* (groups of nodes with identical
symbolic structure) plus a ``node -> template`` map.  Here we need three:
``init`` (start constraint), ``interior`` (the repeated middle nodes), and
``terminal`` (goal constraint, no control / no dynamics).  The interior
template's dynamics/Jacobian/Hessian are generated **once** and reused at every
middle node.

Run:
    conda activate nvblox_env
    export LD_LIBRARY_PATH=/home/lander/install_dir/lib:$LD_LIBRARY_PATH
    python examples/ocp_direct_minimal.py
"""

import numpy as np
import casadi as ca

from fatropgen import OcpProxy, make_stage, generate_solver, DirectOcpSolver

K = 20        # number of nodes
DT = 0.1      # time step


def dynamics(x, u):
    """f(x, u) -> x_{k+1}: integrate the double integrator one step."""
    return ca.vertcat(x[0] + DT * x[1],     # p+ = p + dt*v
                      x[1] + DT * u[0])      # v+ = v + dt*a


def build_proxy() -> OcpProxy:
    # Each template owns its own (x, u) symbols. There are no parameters here,
    # so the shared parameter list is empty ([]).
    # --- interior node: dynamics + cost + acceleration bound ---
    xi, ui = ca.SX.sym("x", 2), ca.SX.sym("u", 1)
    interior = make_stage("interior", xi, ui, params=[],
                          cost=ui[0] ** 2,
                          dynamics=dynamics(xi, ui),
                          g_ineq=ui, lb=[-5.0], ub=[5.0])

    # --- initial node: same, plus the equality "start at rest at p=0" ---
    x0, u0 = ca.SX.sym("x", 2), ca.SX.sym("u", 1)
    init = make_stage("init", x0, u0, params=[],
                      cost=u0[0] ** 2,
                      dynamics=dynamics(x0, u0),
                      g_eq=x0 - ca.DM([0.0, 0.0]),     # x == [0, 0]
                      g_ineq=u0, lb=[-5.0], ub=[5.0])

    # --- terminal node: no control, no dynamics, "end at rest at p=1" ---
    xt, ut = ca.SX.sym("x", 2), ca.SX.sym("u", 0)
    terminal = make_stage("terminal", xt, ut, params=[],
                          cost=ca.SX(0.0),
                          g_eq=xt - ca.DM([1.0, 0.0]))  # x == [1, 0]

    node_to_template = [0] + [1] * (K - 2) + [2]
    return OcpProxy(horizon=K, param_dims=[],
                    templates=[init, interior, terminal],
                    node_to_template=node_to_template)


def main():
    proxy = build_proxy()
    print("Generating + compiling the direct OcpAbstract solver ...")
    result = generate_solver(proxy, "minimal_di", "/tmp/fatropgen/minimal_di",
                             verbose=False)

    solver = DirectOcpSolver(result)
    solver.set_option("tolerance", 1e-8)
    solver.set_option("print_level", 0)
    flag = solver.solve()

    p = np.array([solver.get_x(k)[0] for k in range(K)])
    v = np.array([solver.get_x(k)[1] for k in range(K)])
    a = np.array([solver.get_u(k)[0] for k in range(K - 1)])

    print(f"return flag = {flag} (0 = Success), iterations = {solver.iter_count}")
    print(f"p[0] = {p[0]:.4f}  p[-1] = {p[-1]:.4f}  v[-1] = {v[-1]:.4f}")
    print(f"positions:  {np.round(p, 3)}")

    ok = (flag == 0 and abs(p[0]) < 1e-6 and abs(p[-1] - 1.0) < 1e-6
          and abs(v[-1]) < 1e-6 and np.all(np.abs(a) <= 5.0 + 1e-6))
    print("checks:", "OK" if ok else "FAIL")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
