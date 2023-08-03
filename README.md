<!--
Fatrop - A fast trajectory optimization solver
Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.

This file is part of Fatrop.

Fatrop is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fatrop is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Fatrop.  If not, see <http://www.gnu.org/licenses/>.-->
# FATROP
Fatrop is a constrained nonlinear optimal control problem solver that is fast and achieves a high numerical robustness.

The main features of the solver are:
- high numerical robustness thanks to advanced numerical optimization techniques, inspired by [Ipopt](https://coin-or.github.io/Ipopt/)
- fast by exploiting the optimal control problem structure through a specialized linear solver, based on a [generalized Riccati recursion](https://arxiv.org/abs/2302.14836)
- effective handling of path equality and inequality constraints, without relying on penalty methods
- ability to incorporate exact Lagrangian Hessian information
- ability to be initialized from any, possibly infeasible, solution estimate
- interfaced to [rockit](https://gitlab.kuleuven.be/meco-software/rockit), which is a high-level optimal control problem specification framework, built on top of [CasADi](https://web.casadi.org/)

# Installation instructions
## build and install fatrop
At this moment fatrop is only tested on (Ubuntu) Linux machines.
clone the fatrop repository 

    git clone https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop.git --recursive
    cd fatrop

build and install the fatrop project

    mkdir build
    cd build
    cmake -DBLASFEO_TARGET=X64_AUTOMATIC ..
    make -j

if you want to install fatrop on your system: 
    sudo make install

for non-X64 targets change the blasfeo_target parameter according to the table of https://github.com/giaf/blasfeo
## build and install fatropy

    cd fatropy 
    pip install .

## install rockit with fatropy interface 

    git clone https://gitlab.kuleuven.be/meco-software/rockit.git --recursive
    cd rockit
    pip install .

## examples 

https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_rockit_demo

https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_benchmarks

using fatrop from cpp:

fatrop/executables/RunFatrop.cpp

Developer Lander Vanroye (lander.vanroye@kuleuven.be)

Thanks to all contributors:
- Wilm Decré (python bindings, cmake configuration)
- Ajay Sathya (rockit interface)