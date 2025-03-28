
<p align="center">
    <img src="https://raw.githubusercontent.com/meco-group/fatrop/00a50ff13e6a6a8118951e6e5307f5907dda193c/misc/fatrop_logo.svg" width="400" />
</p>

## FATROP
Fatrop is a constrained nonlinear optimal control problem solver that is fast and achieves a high numerical robustness.

The main features of the solver are:
- high numerical robustness thanks to advanced numerical optimization techniques, inspired by [Ipopt](https://coin-or.github.io/Ipopt/)
- fast by exploiting the optimal control problem structure through a specialized linear solver, based on a [generalized Riccati recursion](https://onlinelibrary.wiley.com/doi/full/10.1002/oca.3064)
- effective handling of path equality and inequality constraints, without relying on penalty methods
- ability to incorporate exact Lagrangian Hessian information
- ability to be initialized from any, possibly infeasible, solution estimate

## Upcoming release

A new version of FATROP is on its way.  
The release is currently in the **beta phase** and available for preview in the [`fatropv1`](https://github.com/meco-group/fatrop/tree/fatropv1) branch.  

## Usage
- fatrop can be used using the "low-level" interface by implementing an OcpAbstract class. See fatrop/ocp/OCPAbstract.hpp
- fatrop is also interfaced with CasADi. A version without blasfeo CPU specialization is distributed through Pypi and conda. A usage example can be found [here](https://github.com/jgillis/fatrop_demo). 


## Citing
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


## Developers

Fatrop is developed by [Lander Vanroye](https://www.lvanroye.github.io) at the [KU Leuven Robotics Research Group](https://www.mech.kuleuven.be/robotics) under supervision of [Wilm Decre](https://www.kuleuven.be/wieiswie/en/person/00052672).

Contributors:
- [Ajay Sathya](https://www.kuleuven.be/wieiswie/en/person/00110259) ([rockit](https://gitlab.kuleuven.be/meco-software/rockit) interface)
- The Fatrop logo was designed by [Louis Callens](https://www.kuleuven.be/wieiswie/en/person/00143705)