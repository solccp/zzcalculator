# Compile the code:
    mkdir /path/to/the/build
    cd /path/to/the/build
    cmake /path/to/the/code
    make
    make test

# Enabling OpenMP parallelism
    -DUSE_OPENMP=ON

# Changing the install path
    -DCMAKE_INSTALL_PREFIX=/path/to/install

#Tested CMake versions:
  - v2.8.2

#Tested Compiler Versions and Platforms
  - gfortran 4.8.2 on Ubuntu 14.04 x86_64
  - Intel Fortran Compiler 12.0 on Ubuntu 10.10 x86_64
  - G95 0.93 on Ubuntu 10.10 x86_64