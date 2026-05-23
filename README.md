## FATROP

Fatrop is a constrained nonlinear optimal control problem solver that is fast and achieves a high numerical robustness.

$$
\begin{align}
\underset{\mathbf{x}_k, \mathbf{u}_k}{\min} \quad & \sum_{k=0}^{K-1} l_k(\mathbf{u}_k, \mathbf{x}_k) \\
\text{subject to} \quad & \mathbf{x}_{k+1} = \mathbf{f}_k(\mathbf{u}_k, \mathbf{x}_k), \quad k = 0, \dots, K-2 \\
& \mathbf{L}_k \leq \mathbf{g}_k(\mathbf{u}_k, \mathbf{x}_k) \leq \mathbf{U}_k \\
& \mathbf{h}_k(\mathbf{u}_k, \mathbf{x}_k) = \mathbf{0}, \quad k = 0, \dots, K-1
\end{align}
$$

The main features of the solver are:

- high numerical robustness thanks to advanced numerical optimization techniques, inspired by [Ipopt](https://coin-or.github.io/Ipopt/)
- fast by exploiting the optimal control problem structure through a specialized linear solver, based on a [generalized Riccati recursion](https://onlinelibrary.wiley.com/doi/full/10.1002/oca.3064)
- high performance linear algebra through integration of [BLASFEO](https://github.com/giaf/blasfeo)
- effective handling of path equality and inequality constraints, without relying on penalty methods
- ability to incorporate exact Lagrangian Hessian information
- ability to be initialized from any, possibly infeasible, solution estimate
- support for optimization over manifolds (Lie groups): primal variables can live on a manifold while the search direction lives in the tangent space, recovered via a user-supplied retraction

## Other problem classes

In addition to the OCP interface above, Fatrop also exposes two general-purpose NLP interfaces that reuse the same interior-point algorithm:

**Dense**: general nonlinear programs, suitable for small-scale dense problems such as inverse kinematics, PnP, and point cloud alignment.

$$
\begin{align}
\underset{\mathbf{x}}{\min} \quad & f(\mathbf{x}) \\
\text{subject to} \quad & \mathbf{g}_\text{eq}(\mathbf{x}) = \mathbf{0} \\
& \mathbf{L} \leq \mathbf{g}_\text{ineq}(\mathbf{x}) \leq \mathbf{U}
\end{align}
$$

**Graph**: block-sparse NLPs whose Hessian sparsity is described by a user-supplied adjacency graph. Suitable for problems such as SLAM and bundle adjustment. Equality constraints are not supported; inequality constraints are block-local (each depending on a single variable block):

$$
\begin{align}
\underset{\mathbf{x}_0, \dots, \mathbf{x}_{N-1}}{\min} \quad & f(\mathbf{x}_0, \dots, \mathbf{x}_{N-1}) \\
\text{subject to} \quad & \mathbf{L}_{k} \leq \mathbf{g}_{k}(\mathbf{x}_k) \leq \mathbf{U}_{k}, \quad k = 0, \dots, N-1
\end{align}
$$

## Getting Started

Refer to the Fatrop [GitHub Pages](https://meco-group.github.io/fatrop/) as a reference to help getting started.

## Usage
- Fatrop can be used using the "low-level" interface by implementing an OcpAbstract class. See [include/fatrop/ocp/ocp_abstract.hpp](include/fatrop/ocp/ocp_abstract.hpp)
- Fatrop is also interfaced with CasADi. A usage example can be found [here](https://github.com/jgillis/fatrop_demo)
- usage examples of both interfaces can be found in the [examples folder](examples)

## Getting Started

Refer to [github pages](https://meco-group.github.io/fatrop/) as a reference to help getting started.

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