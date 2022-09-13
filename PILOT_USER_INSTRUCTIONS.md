# FATROP

Fast Trajectory Optimizer (FATROP) is an efficient and reliable solver for nonlinear optimal control problems with stagewise constraints, aimed at online applications.

This file contains installation and usage instructions for pilot users. Note that Fatrop is still in development, such that these instructions might change without notice.

Instructions are given and tested for Ubuntu 20.04 and 22.04.

# Service desk

If you encounter any problems when installing or using fatrop, you can send an e-mail to the [service desk](mailto:gitlab-incoming+robotgenskill-fatrop-fatrop-5447-issue-@kuleuven.be). This will create an issue, that Lander and I, and the other pilot users can respond to and follow [here](https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop/-/issues/service_desk).

# Installation using binaries

* `pip install --upgrade pip numpy casadi matplotlib`
* `pip install -i https://test.pypi.org/simple/ fatrop`

You can now jump to installing the Fatrop examples.

# Alternative: installation from source

## Install dependencies

### gcc

* `sudo apt install gcc-11`

(Note: this package cannot be found on ubuntu20, probably missing a PPA here. But my gcc version 9.4.0 was sufficient anyway.)

* `sudo apt install gfortran`

### CMake

    sudo apt-get install cmake-curses-gui

### Blasfeo

* go to the directory in which you want to clone Blasfeo, for example: `cd ~/git`
* clone Blasfeo: `git clone https://github.com/giaf/blasfeo.git`
* configure building blasfeo: `cd blasfeo && mkdir -p build && cd build && ccmake ..`
* In ccmake, press `c` to configure initially and set the following: `CMAKE_INSTALL_PREFIX` to `/usr/local` (other options are also possible of course, but this is convenient), and `TARGET` to `X64_AUTOMATIC`, `BUILD_SHARED_LIBS` to `ON`.
* press `c` to configure, `g` to generate
* build blasfeo: `make -j4` (you can change 4 (allowed number of jobs) depending on your CPU)
* install blasfeo: `sudo make install`

### Eigen

* install eigen: `sudo apt install libeigen3-dev`

### Python packages and tools (for fatropy, will however be automatically installed with fatropy as well)

* `pip install pip --upgrade`
* `pip install setuptools --upgrade`
* `pip install Cython`

## Install Fatrop

* go to the directory in which you want to clone fatrop, for example: `cd ~/git`
* clone the fatrop repository: `git clone git@gitlab.kuleuven.be:robotgenskill/fatrop/fatrop.git`
* `cd fatrop`
* switch to the develop branch: `git checkout develop`
* configure building fatrop: `mkdir -p build && cd build && ccmake ..`
* press `c` to configure, `g` to generate
* build fatrop: `make -j4` (you can change 4 (allowed number of jobs) depending on your CPU)
* install fatrop: `sudo make install`

* update links to shared libraries: `sudo ldconfig`

## Install the Fatrop Python packages

* `cd ..` (you should go one directory higher than `build`)
* `cd fatropy`
* optional, but preferred: activate the Python virtual environment you want to use
* `pip install -e .`

## Get and run the Fatrop examples

* go to the directory in which you want to clone the fatrop-examples, for example: `cd ~/git`
* clone the examples: `git clone git@gitlab.kuleuven.be:robotgenskill/fatrop/fatrop-examples.git`
* `cd fatrop-examples`
* run the rocket example: `cd rocket && python Rocket_example.py`
* this solves the example problem with ipopt, generates C-code and a json configuration file that are used by Fatrop, and then solves the problem using Fatrop
* (optionally, you can also solve the problem with fatrop from the terminal: first compile the generated C-code, `gcc -fPIC -march=native -shared -O3 Rocket_example.c -o Rocket_example.so`, and then run Fatrop: `RunFatrop ./Rocket_example.so Rocket_example.json`)

## Additional steps

### Installing Mumps - not needed to use Fatrop, but used for a.o. benchmarking linear solvers

* mumps is not required for fatrop itself, but used for benchmarking.
* install dependencies: `sudo apt install intel-mkl-full libmetis-dev libscotch-dev openmpi-bin libopenmpi-dev`
* go to the directory in which you want to clone mumps, for example: `cd ~/git`
* clone mumps (we use Mumps via CMake): `git clone https://github.com/scivision/mumps`
* configure building mumps: `cd mumps && mkdir -p build && cd build && ccmake ..`
* in `ccmake` set parallel to `OFF`
* if you have a too old version of cmake, install a more recent one with snap: `sudo snap install cmake`
* press `c` to configure, `g` to generate
* build mumps: `make -j4` (you can change 4 (allowed number of jobs) depending on your CPU)
* install mumps: `sudo make install`
