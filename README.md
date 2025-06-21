## FATROP
Fatrop is a constrained nonlinear optimal control problem solver that is fast and achieves a high numerical robustness.

The main features of the solver are:
- high numerical robustness thanks to advanced numerical optimization techniques, inspired by [Ipopt](https://coin-or.github.io/Ipopt/)
- fast by exploiting the optimal control problem structure through a specialized linear solver, based on a [generalized Riccati recursion](https://onlinelibrary.wiley.com/doi/full/10.1002/oca.3064)
- high performance linear algebra through integration of [blasfeo](https://github.com/giaf/blasfeo)
- effective handling of path equality and inequality constraints, without relying on penalty methods
- ability to incorporate exact Lagrangian Hessian information
- ability to be initialized from any, possibly infeasible, solution estimate

## Getting Started

Refer to [github pages](https://meco-group.github.io/fatrop/) as a reference to help getting started.

## Usage
- fatrop can be used using the "low-level" interface by implementing an OcpAbstract class. See [include/fatrop/ocp/ocp_abstract.hpp](include/fatrop/ocp/ocp_abstract.hpp)
- fatrop is also interfaced with CasADi. A usage example can be found [here](https://github.com/jgillis/fatrop_demo)
- usage examples of both interfaces can be found in the [examples folder](examples)

## Citing
To cite Fatrop in your academic work, please use the following reference of the [fatrop paper](https://arxiv.org/abs/2303.16746):

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