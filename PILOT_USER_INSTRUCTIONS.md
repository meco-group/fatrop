# FATROP

Fast Trajectory Optimizer (FATROP) is an efficient and reliable solver for nonlinear optimal control problems with stagewise constraints, aimed at online applications.

This file contains installation and usage instructions for pilot users. Note that Fatrop is still in development, such that these instructions might change without notice.

# Service desk

If you encounter any problems when installing or using fatrop, you can send an e-mail the [service desk](mailto:gitlab-incoming+robotgenskill-fatrop-fatrop-5447-issue-@kuleuven.be). This will create an issue, that Lander and I, and the other pilot users can respond to and follow [here](https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop/-/issues/service_desk).

# Installation

## Installing dependencies

### Blasfeo

* go to the directory in which you want to clone Blasfeo, for example: `cd ~/git`
* clone Blasfeo: `git clone https://github.com/giaf/blasfeo.git`
* configure building blasfeo: `cd blasfeo && mkdir -p build && cd build && ccmake ..`
* In ccmake, set the following: `CMAKE_INSTALL_PREFIX` to `/usr/local` (other options are also possible of course, but this is convenient), and `TARGET` to `X64_AUTOMATIC`
* press `c` to configure, `g` to generate
* build blasfeo: `make -j4` (you can change 4 (allowed number of jobs) depending on your CPU)
* install blasfeo: `sudo make install`

### Mumps - you can skip this if you do not want to benchmark

* mumps is not required for fatrop itself, but used for benchmarking.
* install dependencies: `sudo apt install intel-mkl-full libmetis-dev libscotch-dev openmpi-bin libopenmpi-dev`
* go to the directory in which you want to clone mumps, for example: `cd ~/git`
* clone mumps (we use Mumps via CMake): `git clone https://github.com/scivision/mumps`
* configure building mumps: `cd mumps && mkdir -p build && cd build && ccmake ..`
* in `ccmake` set parallel to `OFF`
* press `c` to configure, `g` to generate
* build mumps: `make -j4` (you can change 4 (allowed number of jobs) depending on your CPU)
* install mumps: `sudo make install`

## Installing Fatrop

* go to the directory in which you want to clone fatrop, for example: `cd ~/git`
* clone the fatrop repository: `git clone git@gitlab.kuleuven.be:robotgenskill/fatrop/fatrop.git`
* configure building fatrop: `cd fatrop && mkdir -p build && cd build && ccmake ..`
* press `c` to configure, `g` to generate
* build fatrop: `make -j4` (you can change 4 (allowed number of jobs) depending on your CPU)
* install fatrop: `sudo make install`

* update links to shared libraries: `sudo ldconfig`
* 

