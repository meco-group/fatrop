FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential cmake git python3 python3-dev python3-pip \
    liblapack-dev gfortran wget pkg-config swig \
    python3-tk \
    && rm -rf /var/lib/apt/lists/*

# Make python3 the default python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Build and install BLASFEO
RUN git clone https://github.com/giaf/blasfeo.git && \
    mkdir blasfeo/build && cd blasfeo/build && \
    cmake .. \
    -DTARGET=X64_AUTOMATIC \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc) && make install

# Build and install FATROP
RUN git clone https://github.com/meco-group/fatrop.git && \
    mkdir fatrop/build && cd fatrop/build && \
    cmake .. \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTS=OFF && \
    make -j$(nproc) && make install

# # # Build and install CasADi
RUN git clone https://github.com/casadi/casadi.git && \
    mkdir casadi/build && cd casadi/build && \
    cmake .. \
    -DWITH_IPOPT=ON -DWITH_BUILD_IPOPT=ON \
    -DWITH_FATROP=ON\
    -DWITH_BUILD_MUMPS=ON -DWITH_BUILD_METIS=ON \
    -DWITH_PYTHON=ON -DWITH_PYTHON3=ON \
    -DPYTHON_PREFIX=$(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())') \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DCMAKE_BUILD_TYPE=Release 
RUN cd casadi/build && make -j$(nproc) && make install

RUN pip install matplotlib numpy scipy