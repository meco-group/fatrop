Installation
============

This page will guide you through the process of installing Fatrop.


Binaries
-------------------------------

1. Using Conda (Recommended):

   Fatrop binaries are available on conda:

   .. code-block:: bash

      conda install libfatrop

   Fatrop is shipped by default with the conda CasADi package.
   To install, simply run:

   .. code-block:: bash

      conda install casadi
   
   .. warning::
      Currently, only Fatrop v0 is available on conda. It is recommended to install the newest version from source for the latest improvements.

   
While it's possible to install Fatrop via PyPI, this method is not recommended as it currently doesn't enable the CPU-specific optimizations.
This means that the performance of the PyPI binaries are significantly lower than the performance of the binaries installed via conda or from source.

Installation from source
-------------------------------

Dependencies
^^^^^^^^^^^^

The only dependency of fatrop is `Blasfeo <https://github.com/giaf/blasfeo>`_.
Blasfeo is a high-perfomance library for linear algebra operations.

1. Clone the blasfeo repository

   .. code-block:: bash

      git clone https://github.com/giaf/blasfeo
      cd blasfeo

2. Create a build directory and navigate to it:

   .. code-block:: bash

      mkdir build
      cd build

3. Configure the project with CMake, specifying the installation prefix:

   .. code-block:: bash

      cmake .. -DTARGET=`blasfeo target` -DCMAKE_INSTALL_PREFIX=/path/to/install/directory
   
   Here, ``blasfeo target`` should be replaced with the target CPU architecture you want to build for.  
   A list of available targets can be found on the
   `BLASFEO GitHub page <https://github.com/giaf/blasfeo?tab=readme-ov-file#supported-computer-architectures>`_.

   Note: You can use `$CONDA_PREFIX` as the installation directory if you're using a Conda environment.

Fatrop
^^^^^^^^^^^^

1. Clone the Fatrop repository:

   .. code-block:: bash

      git clone git@github.com:meco-group/fatrop.git
      cd fatrop

2. Create a build directory and navigate to it:

   .. code-block:: bash

      mkdir build
      cd build

3. Configure the project with CMake, specifying the installation prefix:

   .. code-block:: bash

      cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/directory

   Note: You can use `$CONDA_PREFIX` as the installation directory if you're using a Conda environment:


4. Build the project:

   .. code-block:: bash

      make

5. Install Fatrop:

   .. code-block:: bash

      make install

   This will install Fatrop to the directory specified by `CMAKE_INSTALL_PREFIX`.

   To verify that Fatrop has been installed correctly, you can run one of the example programs provided in the `examples` directory.