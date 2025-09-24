Getting Started
===============

This guide will help you get started with Fatrop, a fast trajectory optimization package.

Different Interfaces
--------------------

Fatrop offers two main interfaces:

1. Low-level Interface
----------------------

The low-level interface is suitable for:

- Power users who want to implement low-overhead function evaluation
- Users who want to interface their own framework

For examples, refer to the C++ examples in the `examples <https://github.com/meco-group/fatrop/tree/main/examples>`_ directory.

2. CasADi Interface
-------------------

The CasADi interface is recommended for most users, offering:

- Easy-to-use symbolic manipulation
- State-of-the-art automatic differentiation
- Access to different solvers beyond fatrop

This interface is recommended for most users.
Refer to the examples and cookbook for usage.

Next Steps
----------

- Explore the `examples  <https://github.com/meco-group/fatrop/tree/main/examples>`_ directory for more detailed examples
- Read the :doc:`fatrop CasADi cookbook <fatrop_casadi_cookbook>`