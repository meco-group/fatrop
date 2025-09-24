Installation
============

This page will guide you through the process of installing Fatrop.

Binaries
--------

Using Conda (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^

Fatrop v0 binaries are available on conda.

.. warning::
   Currently, only Fatrop v0 is available on conda. It is recommended to install the newest version from source for the latest improvements.

Fatrop is shipped by default with the conda CasADi package. To install, simply run:

.. code-block:: bash

   conda install casadi

Alternatively, you can install the `libfatrop` package directly:

.. code-block:: bash

   conda install libfatrop


Using PyPI (Not Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While it's possible to install Fatrop v0 via PyPI, this method is not recommended as it currently doesn't enable CPU-specific optimizations. This means that the performance of the PyPI binaries is significantly lower than the performance of the binaries installed via conda or from source.

Using Docker
------------

A Dockerfile is available in the root directory of the repository. It can be built and run using the following commands, which should be sufficient to try out the Python examples:

.. code-block:: bash

   docker build -t my-fatrop-image .
   docker run -it --rm -v "$PWD":/workspace -w /workspace my-fatrop-image bash

Installation from Source
------------------------

Dependencies
^^^^^^^^^^^^

The only dependency of Fatrop is `Blasfeo <https://github.com/giaf/blasfeo>`_, a high-performance library for linear algebra operations.

To install Blasfeo from source:

1.  Clone the Blasfeo repository:

    .. code-block:: bash

       git clone https://github.com/giaf/blasfeo.git
       cd blasfeo

2.  Create a build directory and navigate into it:

    .. code-block:: bash

       mkdir build
       cd build

3.  Configure the project with CMake, specifying the installation prefix:

    .. code-block:: bash

       cmake .. -DTARGET=<blasfeo target> -DCMAKE_INSTALL_PREFIX=</path/to/install/directory>
    
    Here, ``<blasfeo target>`` should be replaced with the target CPU architecture you want to build for. A list of available targets can be found on the `BLASFEO GitHub page <https://github.com/giaf/blasfeo?tab=readme-ov-file#supported-computer-architectures>`_.

    .. note::
       You can use ``$CONDA_PREFIX`` as the installation directory if you're using a Conda environment.

4.  Build the project:

    .. code-block:: bash

       make -j$(nproc)

5.  Install Blasfeo:

    .. code-block:: bash

       make install

Fatrop
^^^^^^

To install Fatrop from source:

1.  Clone the Fatrop repository:

    .. code-block:: bash

       git clone https://github.com/meco-group/fatrop.git
       cd fatrop

2.  Create a build directory and navigate into it:

    .. code-block:: bash

       mkdir build
       cd build

3.  Configure the project with CMake, specifying the installation prefix:

    .. code-block:: bash

       cmake .. -DCMAKE_INSTALL_PREFIX=</path/to/install/directory>

4.  Build the project:

    .. code-block:: bash

       make -j$(nproc)

5.  Install Fatrop:

    .. code-block:: bash

       make install

    This will install Fatrop to the directory specified by ``CMAKE_INSTALL_PREFIX``.

    To verify that Fatrop has been installed correctly, you can run one of the example programs provided in the ``examples`` directory.

CasADi (Optional)
^^^^^^^^^^^^^^^^^

While not required for Fatrop itself, CasADi is often used alongside Fatrop for modeling and optimization. To install CasADi from source:

1.  Clone the CasADi repository:

    .. code-block:: bash

       git clone https://github.com/casadi/casadi.git
       cd casadi

2.  Create a build directory and navigate into it:

    .. code-block:: bash

       mkdir build
       cd build

3.  Configure the project with CMake, specifying the desired options:

    .. code-block:: bash

       cmake .. \
           -DWITH_IPOPT=ON -DWITH_BUILD_IPOPT=ON \
           -DWITH_BUILD_MUMPS=ON -DWITH_BUILD_METIS=ON \
           -DWITH_FATROP=ON \
           -DWITH_PYTHON=ON -DWITH_PYTHON3=ON \
           -DPYTHON_PREFIX=$(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())') \
           -DCMAKE_INSTALL_PREFIX=</path/to/install/directory> \
           -DCMAKE_BUILD_TYPE=Release

4.  Build the project:

    .. code-block:: bash

       make -j$(nproc)

5.  Install CasADi:

    .. code-block:: bash

       make install
