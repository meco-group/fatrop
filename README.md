<!--
Fatrop - A fast trajectory optimization solver
 Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.

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

<p align="center">
    <img src="https://raw.githubusercontent.com/meco-group/fatrop/00a50ff13e6a6a8118951e6e5307f5907dda193c/misc/fatrop_logo.svg" width="400" />
</p>



# FATROP
Fatrop is a constrained nonlinear optimal control problem solver that is fast and achieves a high numerical robustness.

The main features of the solver are:
- high numerical robustness thanks to advanced numerical optimization techniques, inspired by [Ipopt](https://coin-or.github.io/Ipopt/)
- fast by exploiting the optimal control problem structure through a specialized linear solver, based on a [generalized Riccati recursion](https://onlinelibrary.wiley.com/doi/full/10.1002/oca.3064)
- effective handling of path equality and inequality constraints, without relying on penalty methods
- ability to incorporate exact Lagrangian Hessian information
- ability to be initialized from any, possibly infeasible, solution estimate
- interfaced to [rockit](https://gitlab.kuleuven.be/meco-software/rockit), which is a high-level optimal control problem specification framework, built on top of [CasADi](https://web.casadi.org/)

## Disclaimer

At this moment the easiest way to get specify fatrop problems is by using the [rockit](https://gitlab.kuleuven.be/meco-software/rockit) interface. See [Install rockit with Fatropy interface](#install-rockit-with-fatropy-interface) for installation instructions. The fatrop-rockit-plugin is not very stable yet, and still under development. Apart form the rockit interface, we are working on a ocp specification framework, especially developed for specifying Fatrop problems.

The spectool has several advantages over the rockit interface:
- it supports fatrop problems in all its generality. For example: multi-stage problems are not supported by the fatrop-rockit-interface
- it's possible to specify and solve problems directly from c++, no need for generating the .so and .json seperately from python
- no need for function compilation, function evaluation can be performed in the casadi virtual machine
- some handy features like custom Jacobian and Hessian expressions for more efficient function evaluation then default casadi AD

To compile fatrop with the (beta) spectool make sure to set the CMake flag `-DWITH_SPECTOOL=ON`. Spectool also requires [CasADi](https://github.com/casadi/casadi) to be installed. Currently, we require the casadi/core/function_internal.hpp header, which is installed when CMake flag `-DINSTALL_INTERNAL_HEADERS=ON` is set when installing CasADi.
<!-- Release is expected end of August 2023. -->

# Installation instructions
At this moment Fatrop is mainly tested on Ubuntu Linux machines. There are two installation types: 
- [FATROP](#fatrop)
  - [Disclaimer](#disclaimer)
- [Installation instructions](#installation-instructions)
  - [Build and install Fatrop and Fatropy](#build-and-install-fatrop-and-fatropy)
  - [Build and install Fatrop only](#build-and-install-fatrop-only)
  - [Install rockit with Fatropy interface](#install-rockit-with-fatropy-interface)
  - [Examples](#examples)
- [Citing](#citing)
- [Developers](#developers)

## Build and install Fatrop and Fatropy
(Recursively) clone this repository

    git clone https://github.com/meco-group/fatrop.git --recursive
    cd fatrop


Set the CMake flags, change the BLASFEO target to your system architecture (see table of https://github.com/giaf/blasfeo)

    export CMAKE_ARGS="-DBLASFEO_TARGET=X64_AUTOMATIC -DENABLE_MULTITHREADING=OFF"

Build and install the Fatropy project

    cd fatropy 
    pip install .

Trouble shoot: 
- make sure you're using the newest pip version (pip install --upgrade pip setuptools)

## Build and install Fatrop only
(Recursively) clone this repository

    git clone https://github.com/meco-group/fatrop.git --recursive
    cd fatrop

Build and install the Fatrop project

    mkdir build
    cd build
    cmake -DBLASFEO_TARGET=X64_AUTOMATIC ..
    make -j

If you want to install Fatrop on your system

    sudo make install

For non-x64 targets change the BLASFEO_target parameter according to the table of https://github.com/giaf/blasfeo

## Install rockit with Fatropy interface 

    git clone https://gitlab.kuleuven.be/meco-software/rockit.git
    git clone https://gitlab.kuleuven.be/u0110259/rockit_fatrop_plugin.git ./rockit/rockit/external/fatrop
    cd rockit
    pip install .

## Examples 
https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_rockit_demo

https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_benchmarks

Using Fatrop from C++: check file `fatrop/executables/RunFatrop.cpp`

# Citing
To cite Fatrop in your academic work, please use the following reference:

```
@inproceedings{vanroye2023fatrop,
  title={Fatrop: A fast constrained optimal control problem solver for robot trajectory optimization and control},
  author={Vanroye, Lander and Sathya, Ajay and De Schutter, Joris and Decr{\'e}, Wilm},
  booktitle={2023 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  pages={10036--10043},
  year={2023},
  organization={IEEE}
}
```



# Developers

Fatrop is developed by [Lander Vanroye](https://www.kuleuven.be/wieiswie/en/person/00116913) at the [KU Leuven Robotics Research Group](https://www.mech.kuleuven.be/robotics) under supervision of [Wilm Decre](https://www.kuleuven.be/wieiswie/en/person/00052672).

Contributors:
- [Ajay Sathya](https://www.kuleuven.be/wieiswie/en/person/00110259) ([rockit](https://gitlab.kuleuven.be/meco-software/rockit) interface)
- The Fatrop logo was designed by [Louis Callens](https://www.kuleuven.be/wieiswie/en/person/00143705)