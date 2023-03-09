# FATROP

Fast Trajectory Optimizer (FATROP) is an efficient and reliable solver for nonlinear optimal control problems with stagewise constraints, aimed at online applications.
# Installation instructions
## build and install fatrop
At this moment fatrop is only tested on linux machines.
clone the fatrop repository 

    git clone https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop.git
    cd fatrop

(optionally) checkout the branch with the most recent fatrop version:

    git checkout develop

load the blasfeo submodule

    git submodule update --recursive --init
build and install the fatrop project

    mkdir build
    cd build
    cmake -DBLASFEO_TARGET=X64_AUTOMATIC ..
    make -j
if you want to install fatrop on your system: 
    sudo make install

for non-X64 targets change the blasfeo_target parameter according to the table of https://github.com/giaf/blasfeo
## build and install fatropy
fatropy is the name of the python bindings for fatrop 
install python-dev tools (relace X.X with your python version)

    sudo apt-get install pythonX.X-dev

navigate to the fatropy source folder

    cd <fatrop_source_dir>/fatropy

tell fatropy where to find the libfatrop.so shared library file: uncomment the following line in the setup.py file and replace INSTALLATION FOLDER with the absolute path to the libfatrop.so file. 

either <fatrop_source_dir>/build/fatrop or your installation directory (for example /usr/local/lib)

        # runtime_library_dirs=["INSTALLATION FOLDER"],

alternatively you could set the LD_LIBRARY_PATH to include the folder with the libfatrop.so shared library file

build and install in your python environment

    pip install -e .

## install custom rockit version
Rockit is the framework that is used for specifying optimal control problems in python. 
We had to make a few modifications for the acados interface such that problems are transcribed in the same way for fatrop and acados.
We are currently cleaning up the code and creating a pull request to the orignal acados repository so that this repository can be used to the in the future to run this benchmark.
For now, please use our custom rockit version:

    git clone https://gitlab.kuleuven.be/robotgenskill/fatrop/rockit.git
    cd rockit
    pip install -e .

also install the rockit-fatrop interface

    cd <rockit_source_dir>/rockit/external
    git submodule add https://gitlab.kuleuven.be/u0110259/rockit_fatrop_plugin.git fatrop

## examples 

https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_rockit_demo

https://gitlab.kuleuven.be/robotgenskill/fatrop/fatrop_benchmarks



Contact
Developer Lander Vanroye (lander.vanroye@kuleuven.be)

Thanks to all contributors:
- Wilm Decré for the cmake configuration, pip executable, pilot user instructions, fatropy python interface, ...
- Ajay Sathya for the rockit interface