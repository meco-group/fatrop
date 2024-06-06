### Compilation instructions for fatrop python bindings with spectool

Install dependencies:

        sudo apt install git wget gfortran liblapack-dev pkg-config cmake gcc swig g++ pybind11-dev

(recommended) create and activate a virtual environment:

        python3 -m venv your_env_name
        source your_env_name/bin/activate

Update pip and setuptools:

        pip install --upgrade pip setuptools

find the python site packages directory:

        export PYTHON_SITE_PACKAGES=$(python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")

Install casadi with python bindings from source (note that this will install casadi on your system, consider using a cmake prefix to install it in a specific directory):

    git clone https://github.com/casadi/casadi.git 
    mkdir casadi/build && cd casadi/build
    cmake \
        -DWITH_IPOPT=OFF -DWITH_BUILD_IPOPT=OFF\
        -DWITH_BUILD_MUMPS=OFF -DWITH_BUILD_METIS=OFF\
        -DWITH_OPENMP=OFF -DWITH_THREAD=OFF -DWITH_PYTHON=ON -DWITH_PYTHON3=ON -DPYTHON_PREFIX=$PYTHON_SITE_PACKAGES\
        -DINSTALL_INTERNAL_HEADERS=ON ..
    make -j && sudo make install

Install fatropy with spectool from source:

    git clone https://github.com/meco-group/fatrop.git --recursive
    export CMAKE_ARGS="-DBLASFEO_TARGET=X64_AUTOMATIC -DENABLE_MULTITHREADING=OFF -DWITH_SPECTOOL=ON" && cd fatropy && python -m pip install .
