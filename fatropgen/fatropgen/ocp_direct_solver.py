"""ctypes loader for a solver ``.so`` produced by :mod:`ocp_codegen`.

Wraps the generated C ABI (``ocp_*``) into a small Python object that sets
parameters / options / initial guess, solves, and reads back the per-node primal
solution + iteration count + solve time.
"""

from __future__ import annotations

import ctypes
import os
from typing import List, Optional, Sequence

import numpy as np

from .ocp_codegen import CodegenResult


class DirectOcpSolver:
    def __init__(self, result: CodegenResult):
        self.result = result
        self.proxy = result.proxy
        self._lib = ctypes.CDLL(result.so_path, mode=ctypes.RTLD_GLOBAL)
        self._decl()
        self._h = self._lib.ocp_create()
        if not self._h:
            raise RuntimeError("ocp_create returned null")

    def _decl(self):
        L = self._lib
        c = ctypes
        dptr = c.POINTER(c.c_double)
        L.ocp_create.restype = c.c_void_p
        L.ocp_set_param.argtypes = [c.c_void_p, c.c_int, dptr, c.c_int]
        L.ocp_set_initial.argtypes = [c.c_void_p, dptr, dptr]
        L.ocp_clear_initial.argtypes = [c.c_void_p]
        L.ocp_set_opt_d.argtypes = [c.c_void_p, c.c_char_p, c.c_double]
        L.ocp_set_opt_i.argtypes = [c.c_void_p, c.c_char_p, c.c_longlong]
        L.ocp_set_opt_b.argtypes = [c.c_void_p, c.c_char_p, c.c_int]
        L.ocp_solve.argtypes = [c.c_void_p]; L.ocp_solve.restype = c.c_int
        L.ocp_iter_count.argtypes = [c.c_void_p]; L.ocp_iter_count.restype = c.c_int
        L.ocp_solve_time.argtypes = [c.c_void_p]; L.ocp_solve_time.restype = c.c_double
        L.ocp_horizon.argtypes = [c.c_void_p]; L.ocp_horizon.restype = c.c_int
        L.ocp_nx.argtypes = [c.c_void_p, c.c_int]; L.ocp_nx.restype = c.c_int
        L.ocp_nu.argtypes = [c.c_void_p, c.c_int]; L.ocp_nu.restype = c.c_int
        L.ocp_get_x.argtypes = [c.c_void_p, c.c_int, dptr]
        L.ocp_get_u.argtypes = [c.c_void_p, c.c_int, dptr]
        L.ocp_destroy.argtypes = [c.c_void_p]

    # --- configuration ---
    def set_param(self, idx: int, value):
        v = np.ascontiguousarray(np.asarray(value, dtype=np.float64).reshape(-1, order="F"))
        self._lib.ocp_set_param(self._h, idx, v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), v.size)

    def set_option(self, name: str, value):
        nm = name.encode()
        if isinstance(value, bool):
            self._lib.ocp_set_opt_b(self._h, nm, 1 if value else 0)
        elif isinstance(value, int):
            self._lib.ocp_set_opt_i(self._h, nm, value)
        else:
            self._lib.ocp_set_opt_d(self._h, nm, float(value))

    def set_initial(self, x_per_node: Sequence[np.ndarray], u_per_node: Sequence[np.ndarray]):
        xf = np.ascontiguousarray(np.concatenate([np.asarray(x, float).reshape(-1) for x in x_per_node]))
        uf = np.ascontiguousarray(np.concatenate([np.asarray(u, float).reshape(-1) for u in u_per_node])
                                  if any(len(u) for u in u_per_node) else np.zeros(0))
        dptr = ctypes.POINTER(ctypes.c_double)
        self._lib.ocp_set_initial(self._h, xf.ctypes.data_as(dptr), uf.ctypes.data_as(dptr))

    def clear_initial(self):
        self._lib.ocp_clear_initial(self._h)

    # --- solve / readback ---
    @property
    def horizon(self) -> int:
        return self._lib.ocp_horizon(self._h)

    def nx(self, k: int) -> int:
        return self._lib.ocp_nx(self._h, k)

    def nu(self, k: int) -> int:
        return self._lib.ocp_nu(self._h, k)

    def solve(self) -> int:
        return self._lib.ocp_solve(self._h)

    @property
    def iter_count(self) -> int:
        return self._lib.ocp_iter_count(self._h)

    @property
    def solve_time(self) -> float:
        return self._lib.ocp_solve_time(self._h)

    def get_x(self, k: int) -> np.ndarray:
        n = self.nx(k)
        out = np.zeros(n, dtype=np.float64)
        self._lib.ocp_get_x(self._h, k, out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out

    def get_u(self, k: int) -> np.ndarray:
        n = self.nu(k)
        out = np.zeros(n, dtype=np.float64)
        if n:
            self._lib.ocp_get_u(self._h, k, out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return out

    def x_solution(self) -> List[np.ndarray]:
        return [self.get_x(k) for k in range(self.horizon)]

    def u_solution(self) -> List[np.ndarray]:
        return [self.get_u(k) for k in range(self.horizon)]

    def __del__(self):
        try:
            if getattr(self, "_h", None):
                self._lib.ocp_destroy(self._h)
                self._h = None
        except Exception:
            pass
