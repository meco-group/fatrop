# CasADi Code Generation Example for FATROP

This folder demonstrates code generation for FATROP using CasADi.

## Files

- `generate_casadi.py`: Python script that generates C code (`casadi_generated.c` and `casadi_generated.h`) for a specific optimization problem.
- `ocp_impl_example_codegen.cpp`: C++ example showing how to use the generated C code with FATROP.

## Process

1. `generate_casadi.py` defines an optimal control problem (OCP) using CasADi.
2. The OCP is encapsulated in a CasADi function.
3. CasADi generates C code for this function, including the OCP and its derivatives, along with the underlying solver.
4. `ocp_impl_example_codegen.cpp` interfaces with FATROP using the generated functions.


## Usage

1. Run `python generate_casadi.py` to generate C code.
2. Compile the C++ example with FATROP and the generated code.
3. Run the executable to solve the optimization problem.

For CasADi to be able to detect the problem structure variables and constraints have to be defined in a specific order. More information and examples can be found here: https://github.com/jgillis/fatrop_demo. 

For more information about Casadi codegen usage refer to https://github.com/casadi/casadi/blob/main/docs/examples/cplusplus/codegen_usage.cpp.
