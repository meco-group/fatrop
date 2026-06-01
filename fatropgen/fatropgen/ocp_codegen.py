"""Turn an :class:`~fatropgen.ocp_proxy.OcpProxy` into a
compiled Fatrop solver that implements ``fatrop::Nlp<OcpType>`` directly
(the higher-level interface above ``OcpAbstract``).

Pipeline
--------
1. For each node template build a small set of CasADi ``Function``s in the I/O
   convention ``ocp_direct_runtime.hpp`` expects:

       val_<tid>   : (x, u, params...) -> (b?, g?, gineq?)       [one CFun call]
       jac_<tid>   : (x, u, params...) -> (babt?, ggt?, ggt_ineq?) [one CFun call]
       obj_<tid>   : (x, u, params...) -> (L, rq)
       rsqrqt_<tid>: (x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, params...) -> H
       x0_<tid>/u0_<tid> : (params...) -> initial guess

   The grouped ``val`` and ``jac`` functions are the key shape change versus the
   previous (OcpAbstract) codegen: every constraint Jacobian of a stage now
   comes out of ONE CasADi call instead of three (BAbt + Ggt + Ggt_ineq), which
   matches the natural granularity of ``fatrop::Nlp::eval_constr_jac``.
   Optional ``deriv_overrides["jac"]`` / ``["val"]`` / ``["obj"]`` on a template
   replace the AD-derived versions verbatim.

2. Add every function **once** to a single ``ca.CodeGenerator`` → one C file
   (CasADi dedups shared inner helpers).  A grouped multi-output function makes
   CSE across e.g. ``b`` and ``babt`` (which both use ``qdd = ABA(q,qd,tau)``)
   actually take effect.

3. Render ``<name>_ocp.cpp`` (the ``fatrop::Nlp<OcpType>`` subclass
   instantiation + a tiny C ABI) and compile ``<name>_casadi.c`` +
   ``<name>_ocp.cpp`` into one ``.so``, linking Fatrop/blasfeo (+ any external
   plugin libs).

The generated ``.so`` is loaded by :class:`DirectOcpSolver` (ctypes).

IMPORTANT build note: compile against ``$FATROP_DIR`` headers that match the
runtime ``libfatrop.so`` — mixing the install tree headers with a different
build of the lib causes heap corruption.  Defaults to ``/usr/local`` (matches
the docker container's ``CMAKE_INSTALL_PREFIX``); override via the
``FATROP_DIR`` env var on host installs that use a different prefix
(e.g. ``/home/lander/install_dir``).
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
_KINDS = ["b", "g", "gineq", "jac", "obj", "rsqrqt", "x0", "u0"]
_ENUM_OF = {"b": "B", "g": "G", "gineq": "GINEQ", "jac": "JAC", "obj": "OBJ",
            "rsqrqt": "RSQRQT", "x0": "X0", "u0": "U0"}

_DEFAULT_FATROP_DIR = os.environ.get("FATROP_DIR", "/usr/local")
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


@dataclass
class _TemplateBuild:
    """All the per-template artefacts needed by the renderer."""
    fns: Dict[str, ca.Function]          # name -> Function
    meta: dict                            # dims + lb/ub + i/o indices
    jac_outputs: List[str]                # list of "babt"/"ggt"/"ggt_ineq" present in jac


def _build_template(t: StageProxy, tid: int) -> _TemplateBuild:
    """Build the CasADi Functions for one template."""
    x, u, P = t.x, t.u, list(t.params)
    z = _zcat(u, x)
    nz = z.numel()
    sfx = f"_t{tid}"
    opts = {"cse": True}
    sym = ca.SX.sym if isinstance(x, ca.SX) else ca.MX.sym

    fns: Dict[str, ca.Function] = {}

    def F(name, ins, outs):
        return _try_expand(ca.Function(name, ins, [ca.densify(o) for o in outs], opts))

    # ------------------------------------------------------------------ values
    # The constraint *value* blocks (b, g, gineq) are emitted as separate
    # single-output Functions.  Combining them into one multi-output Function
    # would let CSE share work between blocks that route through the same
    # opaque big-function call (g_all), but CasADi's .expand() trips up on
    # multi-output Functions whose outputs both reference the same opaque
    # CallSX with a nested ca.external (the nvblox SDF case): construction
    # of the forward sensitivity function asserts ``Duplicate``.
    # The Jacobian function (jac) DOES group all stage Jacobians via the
    # `J_all` kernel which is constructed differently (per-template AD lives
    # at the wrapper level, not inside the opaque kernel), so it doesn't hit
    # the same path.
    if t.has_dynamics:
        fns["b"] = F("b" + sfx, [x, u, *P], [t.dynamics])
    if t.ng_eq > 0:
        fns["g"] = F("g" + sfx, [x, u, *P], [t.g_eq])
    if t.ng_ineq > 0:
        fns["gineq"] = F("gineq" + sfx, [x, u, *P], [t.g_ineq])

    # ------------------------------------------------------------------ jac
    # ``jac_outputs`` only describes which blocks are present (so the
    # generated cpp can wire up output indices); the actual ``ca.jacobian``
    # calls are LAZY — only built when no override is provided.  This
    # matters when ``t.g_eq``/``t.g_ineq`` route through an opaque
    # ``g_all`` kernel that contains an ``ca.external`` (e.g. the nvblox
    # SDF): forward-AD through such a kernel asserts ``Duplicate`` during
    # SX-init.  The multi-stage solver always provides ``jac`` precisely
    # to bypass this AD path.
    jac_outputs = []
    if t.has_dynamics:
        jac_outputs.append("babt")
    if t.ng_eq > 0:
        jac_outputs.append("ggt")
    if t.ng_ineq > 0:
        jac_outputs.append("ggt_ineq")
    jac_override = t.deriv_overrides.get("jac")
    if jac_override is not None:
        fns["jac"] = jac_override
        if jac_override.n_out() != len(jac_outputs):
            raise ValueError(
                f"template {t.name!r}: jac override has {jac_override.n_out()} outputs, "
                f"expected {len(jac_outputs)} (babt?+ggt?+ggt_ineq?)")
    elif jac_outputs:
        jac_exprs = []
        if t.has_dynamics:
            jac_exprs.append(ca.jacobian(t.dynamics, z).T)
        if t.ng_eq > 0:
            jac_exprs.append(ca.jacobian(t.g_eq, z).T)
        if t.ng_ineq > 0:
            jac_exprs.append(ca.jacobian(t.g_ineq, z).T)
        fns["jac"] = F("jac" + sfx, [x, u, *P], jac_exprs)

    # ----------------------------------------------------------- cost & grad
    obj_override = t.deriv_overrides.get("obj")
    if obj_override is not None:
        fns["obj"] = obj_override
        if obj_override.n_out() != 2:
            raise ValueError(
                f"template {t.name!r}: obj override must have 2 outputs (L, rq)")
    else:
        # The cost rarely shares structure with the constraint kernels, so
        # the (L, rq) multi-output Function is safe to .expand().
        fns["obj"] = F("obj" + sfx, [x, u, *P], [t.cost, ca.jacobian(t.cost, z).T])

    # ------------------------------------------------------------ rsqrqt
    lam_dyn = sym("lam_dyn" + sfx, t.nx_next)
    lam_eq = sym("lam_eq" + sfx, t.ng_eq)
    lam_ineq = sym("lam_ineq" + sfx, t.ng_ineq)
    obj_scale = sym("obj_scale" + sfx, 1)
    rsqrqt_override = t.deriv_overrides.get("rsqrqt")
    if rsqrqt_override is not None:
        fns["rsqrqt"] = rsqrqt_override
    else:
        lagr = obj_scale * t.cost
        if t.has_dynamics:
            lagr = lagr + ca.dot(lam_dyn, t.dynamics)   # -x_{k+1} drops out of the Hessian
        if t.ng_eq > 0:
            lagr = lagr + ca.dot(lam_eq, t.g_eq)
        if t.ng_ineq > 0:
            lagr = lagr + ca.dot(lam_ineq, t.g_ineq)
        H = ca.hessian(lagr, z)[0]
        fns["rsqrqt"] = F(
            "rsqrqt" + sfx, [x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, *P], [H])

    # ------------------------------------------------------------ initial guess
    if t.x0 is not None:
        fns["x0"] = F("x0" + sfx, [*P], [t.x0])
    if t.u0 is not None and t.nu > 0:
        fns["u0"] = F("u0" + sfx, [*P], [t.u0])

    meta = dict(nx=t.nx, nu=t.nu, nx_next=t.nx_next, ng_eq=t.ng_eq,
                ng_ineq=t.ng_ineq, has_dynamics=t.has_dynamics,
                lb=t.lb.tolist(), ub=t.ub.tolist(),
                jac_outputs=jac_outputs)
    return _TemplateBuild(fns=fns, meta=meta, jac_outputs=jac_outputs)


def _render_cpp(name: str, proxy: OcpProxy, builds: List[_TemplateBuild]) -> str:
    """Render the Nlp<OcpType> subclass instantiation + C ABI."""
    L = []
    w = L.append
    w(f'#include "{name}_casadi.h"')
    w('#include "ocp_direct_runtime.hpp"')
    w('#include <algorithm>')
    w('#include <limits>')
    w('using namespace ocp_direct;')
    w('static std::unique_ptr<CFun> MK(const CasadiCFunctions &f){'
      ' return std::make_unique<CFun>(f); }')
    w('extern "C" {')
    w('void* ocp_create() {')
    w('  auto* s = new Solver();')
    w('  s->nlp = std::make_shared<GeneratedNlp>();')
    w('  auto& O = *s->nlp;')
    w(f'  O.templates.resize({len(builds)});')
    for tid, b in enumerate(builds):
        m = b.meta
        w(f'  {{ auto& T = O.templates[{tid}];')
        w(f'    T.nx={m["nx"]}; T.nu={m["nu"]}; T.nx_next={m["nx_next"]};'
          f' T.ng_eq={m["ng_eq"]}; T.ng_ineq={m["ng_ineq"]};'
          f' T.has_dynamics={"true" if m["has_dynamics"] else "false"};')
        if m["lb"]:
            lb = ",".join(_fmt_double(v) for v in m["lb"])
            ub = ",".join(_fmt_double(v) for v in m["ub"])
            w(f'    T.lb={{{lb}}}; T.ub={{{ub}}};')
        # jac output index mapping (which output of the multi-output jac
        # Function is each Jacobian block).
        for idx, blk in enumerate(b.jac_outputs):
            w(f'    T.jac_idx_{blk}={idx};')
        for kind, fn in b.fns.items():
            sym = fn.name()
            w(f'    T.fn[{_ENUM_OF[kind]}]=MK(CASADI_C_FUNCTIONS({sym}));')
        w('  }')
    n2t = ",".join(str(t) for t in proxy.node_to_template)
    w(f'  O.node_to_template = {{{n2t}}};')
    pd = ",".join(f'{{{r},{c}}}' for (r, c) in proxy.param_dims)
    w(f'  O.param_dims = {{{pd}}};')
    w(f'  O.params.resize({len(proxy.param_dims)});')
    for i, (r, c) in enumerate(proxy.param_dims):
        w(f'  O.params[{i}].assign({int(r) * int(c)}, 0.0);')
    # sensible defaults (overridable from Python)
    w('  s->opt_d["tolerance"]=1e-4; s->opt_i["max_iter"]=200; s->opt_i["print_level"]=0;')
    w('  return s;')
    w('}')
    w('void ocp_set_param(void* h,int idx,const double* v,int n){'
      ' auto* s=(Solver*)h; std::copy(v,v+n,s->nlp->params[idx].begin()); }')
    w('void ocp_set_initial(void* h,const double* xf,const double* uf){'
      ' auto* s=(Solver*)h; auto& O=*s->nlp; int K=O.horizon();'
      ' O.x0_cache.resize(K); O.u0_cache.resize(K); int px=0,pu=0;'
      ' for(int k=0;k<K;++k){ int nx=O.tmpl(k).nx, nu=O.tmpl(k).nu;'
      ' O.x0_cache[k].assign(xf+px,xf+px+nx); px+=nx;'
      ' O.u0_cache[k].assign(uf+pu,uf+pu+nu); pu+=nu; } O.init_overridden=true; }')
    w('void ocp_set_opt_d(void* h,const char* n,double v){ ((Solver*)h)->opt_d[n]=v; }')
    w('void ocp_set_opt_i(void* h,const char* n,long long v){ ((Solver*)h)->opt_i[n]=v; }')
    w('void ocp_set_opt_b(void* h,const char* n,int v){ ((Solver*)h)->opt_b[n]=(v!=0); }')
    w('void ocp_clear_initial(void* h){ ((Solver*)h)->nlp->init_overridden=false; }')
    w('int  ocp_solve(void* h){ return ((Solver*)h)->solve(); }')
    w('int  ocp_iter_count(void* h){ return ((Solver*)h)->iter_count; }')
    w('double ocp_solve_time(void* h){ return ((Solver*)h)->solve_time; }')
    w('int  ocp_horizon(void* h){ return ((Solver*)h)->nlp->horizon(); }')
    w('int  ocp_nx(void* h,int k){ return ((Solver*)h)->nlp->tmpl(k).nx; }')
    w('int  ocp_nu(void* h,int k){ return ((Solver*)h)->nlp->tmpl(k).nu; }')
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
    """Generate + compile a direct-Nlp Fatrop solver from ``proxy``.

    extra_link: list of (lib_dir, lib_name) external plugins to link (e.g. nvblox).
    cache: if True (default), skip the (slow) C/C++ compile when the freshly
        generated sources are byte-identical to those behind an existing ``.so``
        (a content hash is stored next to it). Recompiles automatically whenever
        the generated code changes.
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1. build per-template functions
    builds: List[_TemplateBuild] = []
    for tid, t in enumerate(proxy.templates):
        builds.append(_build_template(t, tid))

    # 2. one CodeGenerator, each function added once (dedup by name)
    cg = ca.CodeGenerator(f"{name}_casadi", {"with_header": True, "with_export": True})
    seen = set()
    for b in builds:
        for fn in b.fns.values():
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
    cpp = _render_cpp(name, proxy, builds)
    cpp_path = os.path.join(out_dir, f"{name}_ocp.cpp")
    with open(cpp_path, "w") as fh:
        fh.write(cpp)

    # 4. compile -> one .so (skipped if the sources are unchanged and cache=True)
    c_path = os.path.join(out_dir, f"{name}_casadi.c")
    so_path = os.path.join(out_dir, f"lib{name}.so")
    hash_path = os.path.join(out_dir, f"{name}.buildhash")
    link = list(extra_link or [])
    # Hash the generated C/C++ AND the header-only runtime + c-wrap, so editing
    # any of them invalidates a cached .so.
    hdr = ""
    for h in ("ocp_direct_runtime.hpp", "casadi_c_wrap.hpp"):
        try:
            hdr += open(os.path.join(_CASADI_UTILS, h)).read()
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
