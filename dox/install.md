Installation {#install}
=======================

Installing the dependencies
============================

Building OpenCAP requires the following:

* C++ compiler with full C++17 language support and standard libraries (**Warning: Default Apple Clang on MacOS requires special care**)

* Python3 version >= 3.4

* [CMake](https://cmake.org/):  version >= 3.12

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/): hierarchical data format, version >= 1.10

* [Eigen](http://eigen.tuxfamily.org/dox/): linear algebra library, version >= 3.3

All of these dependencies are available through standard package managers such as 
[Homebrew](https://brew.sh/), [Conda](https://docs.conda.io/en/latest/), and yum/apt-get 
on Linux. Please consult their respective websites for more information, and ensure that you 
have recent enough versions installed on your machine.

The following packages are automatically built at the CMake step (no action required):

* [Numgrid](https://github.com/dftlibs/numgrid): numerical integration library

* [h5pp](https://github.com/DavidAce/h5pp): C++17 wrapper for HDF5

* [pybind11](https://github.com/pybind/pybind11): C++ Python bindings


Basic steps
============

These steps have been tested on MacOS 10.13-10.15, and on Centos7.

First clone the git repo 

    git clone https://github.com/gayverjr/opencap.git
    cd opencap/opencap

Then generate the Makefile using CMake

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

Once the Makefile is generated, build and install the executable.

    make
    make install

__MacOS__ 

For Mac users on **MacOS 10.14 Mojave** or **MacOS 10.15 Catalina**, the Apple Clang provided by 
XCode is missing some C++17 standard library features, which will cause the CMake step to fail. A drop-in
replacement for the missing std::filesystem can be automatically downloaded and installed by 
passing -DH5PP_DOWNLOAD_METHOD=fetch as an argument to CMake:

    cmake -DH5PP_DOWNLOAD_METHOD=fetch -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

For Mac users on **MacOS 10.13 High Sierra** or earlier, the Apple Clang provided by XCode 
**will not work**. We suggest installing the latest version of GCC (currently 10.2) 
from [Homebrew](https://brew.sh/), and then setting proper environment 
variables for CMake:

    CC=gcc-10 CXX=g++-10 cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

__Linux__

For Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler.

Verify installation
====================

To verify your installation, run the tests in the build directory

    make test

To run OpenCAP, specify an input file from the command line

    /path/to/opencap input_file.in 

See the [examples](https://github.com/gayverjr/opencap/tree/master/examples/opencap) directory 
and our <a href="input.html">input </a> page to help get you started.

