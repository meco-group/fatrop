"""Turn an :class:`~fatropgen.ocp_proxy.OcpProxy` into a
compiled Fatrop solver that implements ``fatrop::OcpAbstract`` directly.

Pipeline
--------
1. For each node template build the CasADi ``Function``s in the exact I/O
   convention ``ocp_direct_runtime.hpp`` expects (values + stagewise Jacobians +
   Lagrangian Hessian).  Missing derivatives are AD-derived here with the
   decision vector ordered ``z = [u; x]`` (controls first) and the Lagrangian
   ``obj_scale*L + lam_dyn*f + lam_eq*g_eq + lam_ineq*g_ineq`` (the ``-x_{k+1}``
   term is constant in ``z`` so it is dropped from the Hessian).
2. Add every function **once** to a single ``ca.CodeGenerator`` → one C file
   (CasADi dedups shared inner helpers).  Because all nodes that share a
   template reference one function, a Jacobian/Hessian reused at many time steps
   is generated exactly once.
3. Render ``<name>_ocp.cpp`` (the ``OcpAbstract`` subclass instantiation + a tiny
   C ABI) and compile ``<name>_casadi.c`` + ``<name>_ocp.cpp`` into one ``.so``,
   linking Fatrop/blasfeo (+ any external plugin libs).

The generated ``.so`` is loaded by :class:`DirectOcpSolver` (ctypes).

IMPORTANT build note: compile against ``$FATROP_DIR`` headers that match the
runtime ``libfatrop.so`` — mixing the install tree headers with a different
build of the lib causes heap corruption.  Defaults to ``/home/lander/install_dir``.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import casadi as ca

from .ocp_proxy import OcpProxy, StageProxy

# Kind order must match the `enum Kind` in ocp_direct_runtime.hpp.
_KINDS = ["b", "babt", "g", "ggt", "gineq", "ggt_ineq", "L", "rq", "rsqrqt", "x0", "u0"]
_KIND_ENUM = {k: e for e, k in enumerate(
    ["B", "BABT", "G", "GGT", "GINEQ", "GGT_INEQ", "L", "RQ", "RSQRQT", "X0", "U0"])}
_ENUM_OF = dict(zip(_KINDS, ["B", "BABT", "G", "GGT", "GINEQ", "GGT_INEQ", "L", "RQ",
                             "RSQRQT", "X0", "U0"]))

_DEFAULT_FATROP_DIR = os.environ.get("FATROP_DIR", "/home/lander/install_dir")
_CASADI_UTILS = os.path.join(os.path.dirname(__file__), "casadi_utils")


def _fmt_double(v: float) -> str:
    """Format a double as a valid C++ literal (handles +/-inf)."""
    v = float(v)
    if v == float("inf"):
        return "std::numeric_limits<double>::infinity()"
    if v == float("-inf"):
        return "-std::numeric_limits<double>::infinity()"
    return repr(v)


def _try_expand(fn: ca.Function) -> ca.Function:
    """Expand MX -> SX for speed/dense C; keep the MX function if it cannot
    (e.g. contains an external)."""
    try:
        return fn.expand()
    except Exception:
        return fn


def _zcat(u, x):
    """z = [u; x] (controls first). Handles empty u (terminal nodes)."""
    if u is None or u.numel() == 0:
        return x
    return ca.vertcat(u, x)


def _build_template_functions(t: StageProxy, tid: int) -> Tuple[Dict[str, ca.Function], dict]:
    """Build the CasADi Functions for one template. Returns (functions, meta)."""
    x, u, P = t.x, t.u, list(t.params)
    z = _zcat(u, x)
    nz = z.numel()
    sfx = f"_t{tid}"
    opts = {"cse": True}
    fns: Dict[str, ca.Function] = {}
    # Multiplier / scale symbols must match the template's symbol type (SX or MX)
    # so the Lagrangian can mix them with the value expressions.
    sym = ca.SX.sym if isinstance(x, ca.SX) else ca.MX.sym

    def F(name, ins, outs):
        return _try_expand(ca.Function(name, ins, [ca.densify(o) for o in outs], opts))

    # --- cost & gradient ---
    fns["L"] = F("L" + sfx, [x, u, *P], [t.cost])
    fns["rq"] = F("rq" + sfx, [x, u, *P], [ca.jacobian(t.cost, z).T])

    # --- dynamics (non-terminal only) ---
    if t.has_dynamics:
        f = t.dynamics
        fns["b"] = F("b" + sfx, [x, u, *P], [f])
        babt = t.deriv_overrides.get("babt")
        fns["babt"] = babt if babt is not None else F(
            "babt" + sfx, [x, u, *P], [ca.jacobian(f, z).T])

    # --- equality constraints ---
    if t.ng_eq > 0:
        fns["g"] = F("g" + sfx, [x, u, *P], [t.g_eq])
        ggt = t.deriv_overrides.get("ggt")
        fns["ggt"] = ggt if ggt is not None else F(
            "ggt" + sfx, [x, u, *P], [ca.jacobian(t.g_eq, z).T])

    # --- inequality constraints ---
    if t.ng_ineq > 0:
        fns["gineq"] = F("gineq" + sfx, [x, u, *P], [t.g_ineq])
        ggti = t.deriv_overrides.get("ggt_ineq")
        fns["ggt_ineq"] = ggti if ggti is not None else F(
            "ggt_ineq" + sfx, [x, u, *P], [ca.jacobian(t.g_ineq, z).T])

    # --- Lagrangian Hessian wrt z ---
    lam_dyn = sym("lam_dyn" + sfx, t.nx_next)
    lam_eq = sym("lam_eq" + sfx, t.ng_eq)
    lam_ineq = sym("lam_ineq" + sfx, t.ng_ineq)
    obj_scale = sym("obj_scale" + sfx, 1)
    rsqrqt = t.deriv_overrides.get("rsqrqt")
    if rsqrqt is not None:
        fns["rsqrqt"] = rsqrqt
    else:
        lagr = obj_scale * t.cost
        if t.has_dynamics:
            lagr = lagr + ca.dot(lam_dyn, t.dynamics)   # -x_{k+1} drops out of the Hessian
        if t.ng_eq > 0:
            lagr = lagr + ca.dot(lam_eq, t.g_eq)
        if t.ng_ineq > 0:
            lagr = lagr + ca.dot(lam_ineq, t.g_ineq)
        H = ca.hessian(lagr, z)[0]
        fns["rsqrqt"] = F("rsqrqt" + sfx, [x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, *P], [H])

    # --- initial guess (optional) ---
    if t.x0 is not None:
        fns["x0"] = F("x0" + sfx, [*P], [t.x0])
    if t.u0 is not None and t.nu > 0:
        fns["u0"] = F("u0" + sfx, [*P], [t.u0])

    meta = dict(nx=t.nx, nu=t.nu, nx_next=t.nx_next, ng_eq=t.ng_eq,
                ng_ineq=t.ng_ineq, has_dynamics=t.has_dynamics,
                lb=t.lb.tolist(), ub=t.ub.tolist())
    return fns, meta


def _render_cpp(name: str, proxy: OcpProxy,
                tmpl_fns: List[Dict[str, ca.Function]],
                tmpl_meta: List[dict]) -> str:
    """Render the OcpAbstract subclass instantiation + C ABI."""
    L = []
    w = L.append
    w(f'#include "{name}_casadi.h"')
    w('#include "ocp_direct_runtime.hpp"')
    w('#include <algorithm>')
    w('#include <limits>')
    w('using namespace ocp_direct;')
    w('extern "C" {')
    w('void* ocp_create() {')
    w('  auto* s = new Solver();')
    w('  s->ocp = std::make_shared<GeneratedOcp>();')
    w('  auto& O = *s->ocp;')
    w(f'  O.templates.resize({len(tmpl_meta)});')
    for tid, (fns, m) in enumerate(zip(tmpl_fns, tmpl_meta)):
        w(f'  {{ auto& T = O.templates[{tid}];')
        w(f'    T.nx={m["nx"]}; T.nu={m["nu"]}; T.nx_next={m["nx_next"]};'
          f' T.ng_eq={m["ng_eq"]}; T.ng_ineq={m["ng_ineq"]};'
          f' T.has_dynamics={"true" if m["has_dynamics"] else "false"};')
        if m["lb"]:
            lb = ",".join(_fmt_double(v) for v in m["lb"])
            ub = ",".join(_fmt_double(v) for v in m["ub"])
            w(f'    T.lb={{{lb}}}; T.ub={{{ub}}};')
        for kind, fn in fns.items():
            sym = fn.name()
            w(f'    T.fn[{_ENUM_OF[kind]}]=MK_CFUN({sym});')
        w('  }')
    n2t = ",".join(str(t) for t in proxy.node_to_template)
    w(f'  O.node_to_template = {{{n2t}}};')
    pd = ",".join(f'{{{r},{c}}}' for (r, c) in proxy.param_dims)
    w(f'  O.param_dims = {{{pd}}};')
    w(f'  O.params.resize({len(proxy.param_dims)});')
    for i, (r, c) in enumerate(proxy.param_dims):
        w(f'  O.params[{i}].assign({int(r) * int(c)}, 0.0);')
    w('  O.allocate_caches();   // pre-size x0/u0 caches (keeps solve() alloc-free)')
    # sensible defaults (overridable from Python)
    w('  s->opt_d["tolerance"]=1e-4; s->opt_i["max_iter"]=200; s->opt_i["print_level"]=0;')
    w('  return s;')
    w('}')
    w('void ocp_set_param(void* h,int idx,const double* v,int n){'
      ' auto* s=(Solver*)h; std::copy(v,v+n,s->ocp->params[idx].begin()); }')
    w('void ocp_set_initial(void* h,const double* xf,const double* uf){'
      ' auto* s=(Solver*)h; auto& O=*s->ocp; int K=O.get_horizon_length();'
      ' O.x0_cache.resize(K); O.u0_cache.resize(K); int px=0,pu=0;'
      ' for(int k=0;k<K;++k){ int nx=O.get_nx(k),nu=O.get_nu(k);'
      ' O.x0_cache[k].assign(xf+px,xf+px+nx); px+=nx;'
      ' O.u0_cache[k].assign(uf+pu,uf+pu+nu); pu+=nu; } O.init_overridden=true; }')
    w('void ocp_set_opt_d(void* h,const char* n,double v){ ((Solver*)h)->opt_d[n]=v; }')
    w('void ocp_set_opt_i(void* h,const char* n,long long v){ ((Solver*)h)->opt_i[n]=v; }')
    w('void ocp_set_opt_b(void* h,const char* n,int v){ ((Solver*)h)->opt_b[n]=(v!=0); }')
    w('void ocp_clear_initial(void* h){ ((Solver*)h)->ocp->init_overridden=false; }')
    w('int  ocp_solve(void* h){ return ((Solver*)h)->solve(); }')
    w('int  ocp_iter_count(void* h){ return ((Solver*)h)->iter_count; }')
    w('double ocp_solve_time(void* h){ return ((Solver*)h)->solve_time; }')
    w('int  ocp_horizon(void* h){ return ((Solver*)h)->ocp->get_horizon_length(); }')
    w('int  ocp_nx(void* h,int k){ return ((Solver*)h)->ocp->get_nx(k); }')
    w('int  ocp_nu(void* h,int k){ return ((Solver*)h)->ocp->get_nu(k); }')
    w('void ocp_get_x(void* h,int k,double* out){'
      ' auto& v=((Solver*)h)->x_sol[k]; std::copy(v.begin(),v.end(),out); }')
    w('void ocp_get_u(void* h,int k,double* out){'
      ' auto& v=((Solver*)h)->u_sol[k]; std::copy(v.begin(),v.end(),out); }')
    w('void ocp_destroy(void* h){ delete (Solver*)h; }')
    w('}')
    return "\n".join(L) + "\n"


@dataclass
class CodegenResult:
    name: str
    out_dir: str
    so_path: str
    proxy: OcpProxy


def generate_solver(proxy: OcpProxy, name: str, out_dir: str, *,
                    extra_link: Optional[List[Tuple[str, str]]] = None,
                    fatrop_dir: str = _DEFAULT_FATROP_DIR,
                    cache: bool = True,
                    verbose: bool = True) -> CodegenResult:
    """Generate + compile a direct-OcpAbstract Fatrop solver from ``proxy``.

    extra_link: list of (lib_dir, lib_name) external plugins to link (e.g. nvblox).
    cache: if True (default), skip the (slow) C/C++ compile when the freshly
        generated sources are byte-identical to those behind an existing ``.so``
        (a content hash is stored next to it). Recompiles automatically whenever
        the generated code changes.
    """
    os.makedirs(out_dir, exist_ok=True)
    # 1. build per-template functions
    tmpl_fns, tmpl_meta = [], []
    for tid, t in enumerate(proxy.templates):
        fns, meta = _build_template_functions(t, tid)
        tmpl_fns.append(fns)
        tmpl_meta.append(meta)

    # 2. one CodeGenerator, each function added once (dedup by name)
    cg = ca.CodeGenerator(f"{name}_casadi", {"with_header": True, "with_export": True})
    seen = set()
    for fns in tmpl_fns:
        for fn in fns.values():
            if fn.name() in seen:
                continue
            seen.add(fn.name())
            cg.add(fn)
    cwd = os.getcwd()
    try:
        os.chdir(out_dir)
        cg.generate()
    finally:
        os.chdir(cwd)

    # 3. render + write the .cpp
    cpp = _render_cpp(name, proxy, tmpl_fns, tmpl_meta)
    cpp_path = os.path.join(out_dir, f"{name}_ocp.cpp")
    with open(cpp_path, "w") as fh:
        fh.write(cpp)

    # 4. compile -> one .so (skipped if the sources are unchanged and cache=True)
    c_path = os.path.join(out_dir, f"{name}_casadi.c")
    so_path = os.path.join(out_dir, f"lib{name}.so")
    hash_path = os.path.join(out_dir, f"{name}.buildhash")
    link = list(extra_link or [])
    # Hash the generated C/C++ AND the header-only runtime, so editing it
    # invalidates a cached .so.
    hdr = ""
    try:
        hdr = open(os.path.join(_CASADI_UTILS, "ocp_direct_runtime.hpp")).read()
    except OSError:
        pass
    src_hash = hashlib.sha256(
        (open(c_path).read() + cpp + hdr + repr(link)).encode()).hexdigest()
    if (cache and os.path.exists(so_path) and os.path.exists(hash_path)
            and open(hash_path).read().strip() == src_hash):
        if verbose:
            print(f"  [cache] sources unchanged -> reusing {so_path}")
        return CodegenResult(name=name, out_dir=out_dir, so_path=so_path, proxy=proxy)

    inc = os.path.join(fatrop_dir, "include")
    lib = os.path.join(fatrop_dir, "lib")
    link_flags = [f"-L{lib}", f"-Wl,-rpath,{lib}", "-lfatrop"]
    for ld, ln in link:
        link_flags += [f"-L{ld}", f"-Wl,-rpath,{ld}", f"-l{ln}"]
    c_obj = os.path.join(out_dir, f"{name}_casadi.o")
    cpp_obj = os.path.join(out_dir, f"{name}_ocp.o")
    cmds = [
        ["gcc", "-O3", "-march=native", "-fPIC", "-c", c_path, "-o", c_obj],
        ["g++", "-O3", "-march=native", "-std=c++17", "-fPIC", "-c", cpp_path, "-o", cpp_obj,
         f"-I{inc}", f"-I{os.path.abspath(_CASADI_UTILS)}", f"-I{out_dir}"],
        ["g++", "-shared", "-o", so_path, c_obj, cpp_obj, *link_flags],
    ]
    for cmd in cmds:
        if verbose:
            print("  $", " ".join(cmd))
        subprocess.run(cmd, check=True)
    with open(hash_path, "w") as fh:
        fh.write(src_hash)

    return CodegenResult(name=name, out_dir=out_dir, so_path=so_path, proxy=proxy)
