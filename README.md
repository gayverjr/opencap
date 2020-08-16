<div align="left">
  <img src="https://github.com/gayverjr/opencap/blob/master/images/opencap_logo.png" height="200px">
</div>

# OpenCAP: An open-source program for studying resonances in molecules.                                                               
| **Documentation** | [![][docs-img]][docs-url] [![][rtd-img]][rtd-url] |
| :------ | :------- |
| **Build Status**  | [![][travis-img]][travis-url] [![][cov-img]][cov-url] [![][lgtm-img]][lgtm-url]  |

[docs-img]: https://img.shields.io/badge/doc-latest-blue.svg
[docs-url]: https://gayverjr.github.io/opencap/
[travis-img]: https://travis-ci.com/gayverjr/opencap.svg?branch=master
[travis-url]: https://travis-ci.com/gayverjr/opencap
[lgtm-img]: https://img.shields.io/lgtm/grade/cpp/g/gayverjr/opencap.svg?logo=lgtm&logoWidth=18
[lgtm-url]: https://lgtm.com/projects/g/gayverjr/opencap/context:cpp
[rtd-img]: https://readthedocs.org/projects/gayverjropencap/badge/?version=latest
[rtd-url]: https://gayverjropencap.readthedocs.io/en/latest/?badge=latest
[cov-img]: https://codecov.io/gh/gayverjr/opencap/branch/master/graph/badge.svg
[cov-url]: https://codecov.io/gh/gayverjr/opencap/branch/master

OpenCAP is an open-source application aimed at extending the capabilities of electronic structure packages to describe metastable electronic states. 
We currently support an interface with the [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas) and [Pyscf](https://github.com/pyscf/pyscf) packages to 
compute resonance positions and widths using the complex absorbing potential method (CAP). 

Please see the [examples](https://github.com/gayverjr/opencap/tree/master/examples) directory to help get you started on using the software.

# Documentation
Our documentation is currently under construction. We have separate documentation for the 
command line version (which is aimed more at developers) and the Python API (which is aimed 
more at users).

* Doxygen Documentation for C++ code: https://gayverjr.github.io/opencap/

* Readthedocs Documentation for Python API: https://gayverjropencap.readthedocs.io/en/latest/

# Installation

## PyOpenCAP (Python module) 

### Installing with pip (coming soon...)

    pip install pyopencap

Precompiled Python wheels will soon be available on Pypi for almost all Linux systems, and 
most MacOS systems, for Python versions 3.4 and later. This will be the recommended way of 
installing PyOpenCAP. 

### Compiling from source

Compiling PyOpenCAP from source requires the [dependencies](#Dependencies) listed below. 
See our [documentation](https://gayverjropencap.readthedocs.io/en/latest/) for more details. 

If your operating system/Python environment is not covered by any of our pre-built wheels,
pip will then try to compile from source. You can also clone the repository and install 
yourself with pip:

```
git clone https://github.com/gayverjr/opencap.git

cd opencap

pip install .
```

For Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler.

For Mac users, the Apple Clang provided by XCode will not work due to missing standard 
library features. We suggest installing the latest version of GCC (currently 10.2) 
from [Homebrew](https://brew.sh/), and then setting the following 
environment variables before attempting to pip install:

```
# for GCC 10 installed by brew

export CC=gcc-10

export CXX=g++-10
```

To ensure that the installation was successful, start a Python shell, and type:

    import pyopencap

See the [examples](https://github.com/gayverjr/opencap/tree/master/examples/pyopencap) directory 
and our [tutorial](https://gayverjropencap.readthedocs.io/en/latest/tutorial.html) to help get you started.

## OpenCAP (command line version)

Compiling OpenCAP requires the [dependencies](#Dependencies) listed below. See our 
[documentation](https://gayverjr.github.io/opencap/) for more details. Assuming these 
are all installed in locations visible to CMake, installation can proceed as follows:

First clone the git repo

```

git clone https://github.com/gayverjr/opencap.git

cd opencap/opencap

```

Then generate the Makefile using CMake

```

mkdir build

cd build

cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

```

### MacOS 

For Mac users on MacOS 10.14 Mojave or MacOS 10.15 Catalina, the Apple Clang provided by 
XCode is missing some C++17 standard library features, which will cause CMake to fail. A drop-in
replacement for the missing std::filesystem can be automatically downloaded and installed by 
passing -DH5PP_DOWNLOAD_METHOD=fetch as an argument to CMake:

    cmake -DH5PP_DOWNLOAD_METHOD=fetch -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

For Mac users on MacOS 10.13 High Sierra or earlier, the Apple Clang provided by XCode 
will not work. We suggest installing the latest version of GCC (currently 10.2) 
from [Homebrew](https://brew.sh/), and then setting proper environment 
variables for CMake:

    CC=gcc-10 CXX=g++-10 cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..


### Linux

For Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler.


Once the Makefile is generated, build install the executable.

```
make

make test

make install

```

To verify your installation, run the tests in the build directory

    make test

To run OpenCAP, specify an input file from the command line

    /path/to/opencap input_file.in 

See the [examples](https://github.com/gayverjr/opencap/tree/master/examples/opencap) directory 
and our [documentation](https://gayverjr.github.io/opencap/input.html) to help get you started.


# Dependencies

Building OpenCAP/PyOpenCAP from source requires the following:

* C++ compiler with full C++17 language support and standard libraries (*Warning: Apple Clang on MacOS is not supported*)

* Python3 version >= 3.4

* [CMake](https://cmake.org/):  verison >= 3.12

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/): hierarchical data format, version >= 1.8

* [Eigen](http://eigen.tuxfamily.org/dox/): linear algebra library, version >= 3.3


The following packages are automatically built at the CMake step (no action required):

* [Numgrid](https://github.com/dftlibs/numgrid): numerical integration library

* [h5pp](https://github.com/DavidAce/h5pp): C++17 wrapper for HDF5

* [pybind11](https://github.com/pybind/pybind11): C++ Python bindings
 
# Acknowledgements
This project is funded by the Molecular Sciences Software Institute.

![MolSSI logo](https://github.com/gayverjr/OpenCAP/blob/master/images/molssi_logo.png)
