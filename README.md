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
At this moment Fatrop is mainly tested on Ubuntu Linux machines. There are two installation types: 
- [installation of Fatrop (C++) and Fatropy providing Python bindings (using pybind11)](#build-and-install-fatrop-and-fatropy)
- [installation of Fatrop (C++) only](#build-and-install-fatrop-only)

## Build and install Fatrop and Fatropy
(Recursively) clone this repository

    git clone https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop.git --recursive
    cd fatrop

Install dependencies/tools

    sudo apt install pybind11-dev
    pip install --upgrade pip setuptools

Set the CMake flags, change the BLASFEO target to your system architecture (see table of https://github.com/giaf/blasfeo)

    export CMAKE_FLAGS="-DBLASFEO_TARGET=X64_AUTOMATIC -DENABLE_MULTITHREADING=OFF"

Build and install the fatropy project

    cd fatropy 
    pip install .

## Build and install Fatrop only
(Recursively) clone this repository

    git clone https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop.git --recursive
    cd fatrop

<<<<<<< HEAD
Build and install the Fatrop project
=======
## build and install fatrop (c++ only)
build and install the fatrop project
>>>>>>> a05327fe781d9a23b417e3ebc414a707c452f4b9

    mkdir build
    cd build
    cmake -DBLASFEO_TARGET=X64_AUTOMATIC ..
    make -j

If you want to install Fatrop on your system

    sudo make install

For non-x64 targets change the BLASFEO_target parameter according to the table of https://github.com/giaf/blasfeo

## Install rockit with fatropy interface 

    git clone https://gitlab.kuleuven.be/meco-software/rockit.git
    git clone https://gitlab.kuleuven.be/u0110259/rockit_fatrop_plugin.git ./rockit/rockit/external/fatrop --recursive
    cd rockit
    pip install .

## Examples 
https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_rockit_demo

https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_benchmarks

Using Fatrop from C++: check file `fatrop/executables/RunFatrop.cpp`

# Developers

Fatrop is developed by [Lander Vanroye](https://www.kuleuven.be/wieiswie/en/person/00116913) at the [KU Leuven Robotics Research Group](https://www.mech.kuleuven.be/robotics) under supervision of [Wilm Decre](https://www.kuleuven.be/wieiswie/en/person/00052672).

Contributors:
- [Ajay Sathya](https://www.kuleuven.be/wieiswie/en/person/00110259) ([rockit](https://gitlab.kuleuven.be/meco-software/rockit) interface)