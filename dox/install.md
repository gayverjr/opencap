Installation {#install}
=======================

Installing the dependencies
============================

Building OpenCAP requires the following:

* C++ compiler with full C++17 language support and standard libraries

* Python3 interpreter and development libraries: version >= 3.6

* [CMake](https://cmake.org/):  version >= 3.12

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/): hierarchical data format, version >= 1.10

* [Eigen](http://eigen.tuxfamily.org/dox/): linear algebra library, version >= 3.3

All of these dependencies are available through standard package managers such as 
[Homebrew](https://brew.sh/), [Conda](https://docs.conda.io/en/latest/), and yum/apt-get 
on Linux. 

The following packages are automatically built at the CMake step (no action required):

* [Numgrid](https://github.com/dftlibs/numgrid): numerical integration library

* [h5pp](https://github.com/DavidAce/h5pp): C++17 wrapper for HDF5

* [pybind11](https://github.com/pybind/pybind11): C++ Python bindings


Basic steps
============

These steps have been tested on MacOS 12.3, and on Centos7.

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

__MacOS/Linux__ 

For Mac/Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler. 

For Linux users, depending on your distribution, you may also need to install 
Python3 development libraries e.g. `sudo apt-get install python3.x-dev`.

Verify installation
====================

To verify your installation, run the tests in the build directory

    make test

To run OpenCAP, specify an input file from the command line

    /path/to/opencap input_file.in 

See the [examples](https://github.com/gayverjr/opencap/tree/main/examples/opencap) directory 
and our <a href="input.html">input </a> page to help get you started.

