# fatropgen

A small, self-contained Python package that turns a CasADi-described OCP into a
compiled Fatrop solver implementing `fatrop::OcpAbstract` directly (no `Opti` /
`nlpsol` layer).

## Install

```bash
pip install -e .
```

## Native dependencies

The codegen step emits C + C++ and links against **Fatrop + BLASFEO** at build
time, and `dlopen`s them at run time. Set `FATROP_DIR` to the install prefix
that contains `include/fatrop` and `lib/libfatrop.so`:

```bash
export FATROP_DIR=/path/to/fatrop/install
export LD_LIBRARY_PATH=$FATROP_DIR/lib:$LD_LIBRARY_PATH
```

CasADi (with the Fatrop interface compiled in) must also be importable from
Python.

## Run the minimal example

```bash
export LD_LIBRARY_PATH=$FATROP_DIR/lib:$LD_LIBRARY_PATH
python examples/ocp_direct_minimal.py
```

It builds a 1-D double-integrator point-to-point problem, generates + compiles
the solver `.so` under `/tmp/fatropgen/minimal_di/`, and prints the recovered
position trajectory.

## API

```python
from fatropgen import OcpProxy, make_stage, generate_solver, DirectOcpSolver
```

- `make_stage(name, x, u, params, cost, *, dynamics=..., g_eq=..., g_ineq=..., lb=..., ub=...)`
  — symbolic description of one node template.
- `OcpProxy(horizon, param_dims, templates, node_to_template)` — the full OCP.
- `generate_solver(proxy, name, out_dir)` — emits + compiles the solver `.so`.
- `DirectOcpSolver(result)` — `ctypes` wrapper around the generated `.so`;
  `set_option`, `set_param`, `set_initial`, `solve`, `get_x`, `get_u`.

See `examples/ocp_direct_minimal.py` for the smallest end-to-end usage.

## Layout

```
fatropgen/
├── pyproject.toml
├── README.md
├── fatropgen/                # the importable Python package
│   ├── __init__.py
│   ├── ocp_proxy.py          # OcpProxy / StageProxy / make_stage
│   ├── ocp_codegen.py        # generate_solver: CasADi C codegen + g++ link
│   ├── ocp_direct_solver.py  # DirectOcpSolver: ctypes wrapper
│   └── casadi_utils/         # C++ headers, installed as package data
│       └── ocp_direct_runtime.hpp
└── examples/
    └── ocp_direct_minimal.py
```
