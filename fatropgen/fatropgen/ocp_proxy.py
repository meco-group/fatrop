"""A *super-simple* Python+CasADi description of a Fatrop OCP.

This module defines a high-level "proxy" for an optimal-control problem in the
exact shape Fatrop's :class:`fatrop::OcpAbstract` expects, but expressed purely
with CasADi symbols / expressions.  The companion :mod:`ocp_codegen` turns a
proxy into a C++ ``OcpAbstract`` subclass that calls CasADi-generated C code and
solves with Fatrop v1.0.1 directly (no CasADi ``Opti``/``nlpsol`` layer).

Design
------
The OCP is described per **node template** (:class:`StageProxy`) — a group of
nodes that share an identical symbolic structure.  A horizon of ``K`` nodes maps
to a small number of templates via ``node_to_template``; codegen emits one set
of C functions per template, so a Jacobian/Hessian reused at many time steps is
generated exactly once.

Fatrop conventions baked in (see ``fatrop_runtime.hpp`` of the CasADi fatrop
interface, the byte-level reference):

* The per-node decision vector is ``z = [u_k; x_k]`` (**controls first**).
* Dynamics are written as the prediction ``f(x_k, u_k) -> x_{k+1}`` of size
  ``nx_{k+1}`` (which may differ from ``nx_k`` — e.g. a free-time ``T`` that is a
  *control* at node 0 and a *state* afterwards).  Fatrop's defect is
  ``b = -x_{k+1} + f``.
* General constraints are split by Fatrop into equalities ``g_eq == 0`` and
  inequalities ``lb <= g_ineq <= ub`` (constant bounds).

Each :class:`StageProxy` provides the *values* (dynamics ``f``, cost ``L``,
``g_eq``, ``g_ineq``) and the constant ``lb``/``ub``.  Stagewise Jacobians and
the Lagrangian Hessian are derived by CasADi AD in codegen unless explicitly
overridden (see :attr:`StageProxy.deriv_overrides`), which is how a shared
sub-computation (e.g. one ABA call feeding both the dynamics and an acceleration
constraint) is expressed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple, Dict

import numpy as np
import casadi as ca


def _ncol_vec(expr) -> "ca.SX | ca.MX":
    """Return ``expr`` as a single column vector (CasADi ``vertcat`` of a list,
    or the expression itself reshaped to a column)."""
    if expr is None:
        return None
    if isinstance(expr, (list, tuple)):
        if len(expr) == 0:
            return None
        return ca.vvcat(list(expr))
    return ca.reshape(expr, expr.numel(), 1)


@dataclass
class StageProxy:
    """Symbolic description of one node template.

    All expressions are written in terms of the template-local symbols
    :attr:`x` (state, ``nx`` x 1) and :attr:`u` (control, ``nu`` x 1) plus the
    OCP-global parameter symbols :attr:`params` (shared across templates).
    """

    name: str
    x: "ca.SX | ca.MX"          # state symbol, shape (nx, 1)
    u: "ca.SX | ca.MX"          # control symbol, shape (nu, 1); may be empty
    params: List["ca.SX | ca.MX"]  # OCP-global parameter symbols (shared)

    # --- values (expressions) ---
    cost: "ca.SX | ca.MX"                 # scalar stage cost L(x,u,p)
    dynamics: Optional["ca.SX | ca.MX"] = None   # f(x,u,p) -> x_{k+1}; None at terminal
    g_eq: Optional["ca.SX | ca.MX"] = None       # equality value g_eq(x,u,p) == 0
    g_ineq: Optional["ca.SX | ca.MX"] = None     # inequality value, bounded by lb/ub

    # --- constant inequality bounds (numpy arrays sized like g_ineq) ---
    lb: Optional[np.ndarray] = None
    ub: Optional[np.ndarray] = None

    # --- initial guess (expressions in params, or None -> zeros) ---
    x0: Optional["ca.SX | ca.MX"] = None  # (nx, 1)
    u0: Optional["ca.SX | ca.MX"] = None  # (nu, 1)

    # --- optional derivative overrides ---
    # Map from kind in {"val", "jac", "obj", "rsqrqt"} to a pre-built
    # ca.Function with the exact codegen I/O signature (see ocp_codegen).
    # Multi-output kinds (``val``/``jac``) must produce outputs in the order
    # described in ocp_codegen (b? + g? + gineq? for val; babt? + ggt? + ggt_ineq?
    # for jac; the ``?`` block is absent when its dimension is zero).  When an
    # override is present the codegen uses it verbatim instead of AD-deriving.
    deriv_overrides: Dict[str, ca.Function] = field(default_factory=dict)

    def __post_init__(self):
        self.cost = ca.reshape(self.cost, 1, 1)
        self.dynamics = _ncol_vec(self.dynamics)
        self.g_eq = _ncol_vec(self.g_eq)
        self.g_ineq = _ncol_vec(self.g_ineq)
        if self.g_ineq is not None:
            n = self.g_ineq.numel()
            if self.lb is None or self.ub is None:
                raise ValueError(f"template {self.name!r}: g_ineq given without lb/ub")
            self.lb = np.asarray(self.lb, dtype=float).reshape(-1)
            self.ub = np.asarray(self.ub, dtype=float).reshape(-1)
            if self.lb.size != n or self.ub.size != n:
                raise ValueError(
                    f"template {self.name!r}: lb/ub size {self.lb.size}/{self.ub.size} "
                    f"!= g_ineq size {n}")
        else:
            self.lb = np.zeros(0)
            self.ub = np.zeros(0)

    # --- convenient dimensions ---
    @property
    def nx(self) -> int:
        return self.x.numel()

    @property
    def nu(self) -> int:
        return 0 if self.u is None else self.u.numel()

    @property
    def nx_next(self) -> int:
        return 0 if self.dynamics is None else self.dynamics.numel()

    @property
    def ng_eq(self) -> int:
        return 0 if self.g_eq is None else self.g_eq.numel()

    @property
    def ng_ineq(self) -> int:
        return 0 if self.g_ineq is None else self.g_ineq.numel()

    @property
    def has_dynamics(self) -> bool:
        return self.dynamics is not None


@dataclass
class OcpProxy:
    """A complete OCP: a horizon, parameter dimensions, the unique node
    templates, and the node -> template assignment."""

    horizon: int                       # number of nodes K
    param_dims: List[Tuple[int, int]]  # parameter (rows, cols), like get_number_of_parameters()
    templates: List[StageProxy]
    node_to_template: List[int]        # length == horizon

    def __post_init__(self):
        if len(self.node_to_template) != self.horizon:
            raise ValueError(
                f"node_to_template has {len(self.node_to_template)} entries, "
                f"expected horizon={self.horizon}")
        for k, t in enumerate(self.node_to_template):
            if not (0 <= t < len(self.templates)):
                raise ValueError(f"node {k}: template index {t} out of range")
        # The last node must have no dynamics; every other node must have it
        # (Fatrop's shooting structure links node k -> k+1 for k < K-1).
        for k, t in enumerate(self.node_to_template):
            tmpl = self.templates[t]
            if k == self.horizon - 1 and tmpl.has_dynamics:
                raise ValueError(f"terminal node {k} (template {tmpl.name!r}) must not have dynamics")
            if k < self.horizon - 1 and not tmpl.has_dynamics:
                raise ValueError(f"non-terminal node {k} (template {tmpl.name!r}) needs dynamics")
        # Dynamics output size at node k must equal nx at node k+1.
        for k in range(self.horizon - 1):
            t, tn = self.node_to_template[k], self.node_to_template[k + 1]
            got, want = self.templates[t].nx_next, self.templates[tn].nx
            if got != want:
                raise ValueError(
                    f"node {k}: dynamics output dim {got} != node {k+1} nx {want}")

    def template_of(self, k: int) -> StageProxy:
        return self.templates[self.node_to_template[k]]


def make_stage(name, x, u, params, cost, *, dynamics=None, g_eq=None,
               g_ineq=None, lb=None, ub=None, x0=None, u0=None,
               deriv_overrides=None) -> StageProxy:
    """Thin keyword constructor for :class:`StageProxy` (mirrors its fields)."""
    return StageProxy(
        name=name, x=x, u=u, params=list(params), cost=cost, dynamics=dynamics,
        g_eq=g_eq, g_ineq=g_ineq, lb=lb, ub=ub, x0=x0, u0=u0,
        deriv_overrides=dict(deriv_overrides or {}))
