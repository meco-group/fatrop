Installation
============

This page will guide you through the process of installing Fatrop.


Binaries
-------------------------------

If you plan to use Fatrop from the CasADi interface, there are two main methods for installation:

1. Using Conda (Recommended):

   Fatrop binaries are available on conda:

   .. code-block:: bash

      conda install libfatrop


   .. note::
      Currently, only Fatrop v0 is available on conda. It is recommended to install the newest version from source for the latest features and improvements.
   
   Fatrop is shipped by default with the conda CasADi package.
   To install, simply run:

   .. code-block:: bash

      conda install casadi
   


2. Using PyPI:
   
   While it's possible to install Fatrop via PyPI, this method is not recommended as it doesn't enable the CPU-specific optimizations:

Installation from source
-------------------------------

Before installing fatrop, ensure you have the following prerequisites:

BLASFEO can be found at `https://github.com/giaf/blasfeo`.
Alternatively, binaries are available in conda as `libblasfeo`.

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


4. Clone the Fatrop repository:

   .. code-block:: bash

      git clone https://github.com/lvanroye/fatrop.git
      cd fatrop

5. Create a build directory and navigate to it:

   .. code-block:: bash

      mkdir build
      cd build

6. Configure the project with CMake, specifying the installation prefix:

   .. code-block:: bash

      cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/directory

   Note: You can use `$CONDA_PREFIX` as the installation directory if you're using a Conda environment:


7. Build the project:

   .. code-block:: bash

      make

8. Install Fatrop:

   .. code-block:: bash

      make install

   This will install Fatrop to the directory specified by `CMAKE_INSTALL_PREFIX`.

   To verify that Fatrop has been installed correctly, you can run one of the example programs provided in the `examples` directory.
