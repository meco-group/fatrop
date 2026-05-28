"""fatropgen: generate compiled Fatrop OCP solvers from CasADi expressions.

Public API::

    from fatropgen import OcpProxy, StageProxy, make_stage
    from fatropgen import generate_solver, DirectOcpSolver

See ``examples/ocp_direct_minimal.py`` for the smallest end-to-end usage.
"""

from .ocp_proxy import OcpProxy, StageProxy, make_stage
from .ocp_codegen import generate_solver, CodegenResult
from .ocp_direct_solver import DirectOcpSolver

__version__ = "0.1.0"

__all__ = [
    "OcpProxy",
    "StageProxy",
    "make_stage",
    "generate_solver",
    "CodegenResult",
    "DirectOcpSolver",
]
