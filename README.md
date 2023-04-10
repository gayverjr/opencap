<div align="left">
  <img src="https://github.com/gayverjr/opencap/blob/main/images/opencap_logo.png" height="400px">
</div>

# OpenCAP: An open-source program for studying resonances in molecules.                                                               
| **Documentation** | [![][docs-img]][docs-url] [![][rtd-img]][rtd-url] |
| :------ | :------- |
| **Build Status**  | [![][gh-actions-img]][gh-actions-url] [![][cov-img]][cov-url] [![][codeql-img]][codeql-url] |
| **Releases** | [![][pypi-img]][pypi-url] [![][ghrel-img]][ghrel-url] |

[docs-img]: https://img.shields.io/badge/doc-latest-blue.svg
[docs-url]: https://gayverjr.github.io/opencap/
[gh-actions-img]: https://github.com/gayverjr/opencap/workflows/ubuntu/badge.svg
[gh-actions-url]: https://github.com/gayverjr/opencap/actions
[codeql-img]: https://github.com/gayverjr/opencap/actions/workflows/codeql.yml/badge.svg
[codeql-url]: https://github.com/gayverjr/opencap/actions
[rtd-img]: https://readthedocs.org/projects/gayverjropencap/badge/?version=latest
[rtd-url]: https://gayverjropencap.readthedocs.io/en/latest/?badge=latest
[cov-img]: https://codecov.io/gh/gayverjr/opencap/branch/main/graph/badge.svg?token=V9quaajaZ1
[cov-url]: https://codecov.io/gh/gayverjr/opencap/branch/main
[pypi-img]: https://badge.fury.io/py/pyopencap.svg
[pypi-url]: https://pypi.org/project/pyopencap/
[ghrel-img]: https://img.shields.io/github/v/release/gayverjr/opencap.svg
[ghrel-url]: https://github.com/gayverjr/opencap

OpenCAP is an open-source application aimed at extending the capabilities of electronic structure packages to describe metastable electronic states. 
We currently support an interface with the [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas), [Pyscf](https://github.com/pyscf/pyscf), [Q-Chem](https://www.q-chem.com/), [Columbus](https://www.univie.ac.at/columbus/), and [PSI4](https://psicode.org/) packages to 
compute resonance positions and widths using the complex absorbing potential method (CAP). 

Please see the [examples](https://github.com/gayverjr/opencap/tree/main/examples) directory or our [Getting Started Page](https://gayverjropencap.readthedocs.io/en/latest/getting_started.html)
to help get started on using the software. 

For questions or support, please open an issue on GitHub, or contact us directly at gayverjr@bu.edu.

# Documentation
We have separate documentation for the command line version (which is aimed more at developers) and the Python API (which is aimed 
more at users).

* Doxygen Documentation for C++ code: https://gayverjr.github.io/opencap/

* Readthedocs Documentation for Python API: https://gayverjropencap.readthedocs.io/en/latest/

# Installation

## PyOpenCAP (Python module) 

### Install with pip (recommended)

    pip install pyopencap
    # or
    pip3 install pyopencap

Precompiled Python wheels are available on Pypi for almost all Linux systems, and 
most MacOS systems, for Python versions 3.6 and later. 

### Build from source

Compiling PyOpenCAP from source requires first installing all of our [dependencies](#Dependencies). 

For Mac/Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler. 

If your operating system/Python environment is not covered by any of our pre-built wheels,
the command `pip install pyopencap` will download the tarball from Pypi and try to compile from source.
You can also clone the repository and install a local version:

```
git clone https://github.com/gayverjr/opencap.git

cd opencap

pip install .
```

Compiling from source will take several minutes. To monitor your progress, you can run pip 
with the `--verbose` flag.

### Verify installation

To ensure that the installation was successful, return to your home directory, start a Python shell, and type:

    import pyopencap

See the [examples](https://github.com/gayverjr/opencap/tree/main/examples/pyopencap) directory 
and our [tutorial](https://gayverjropencap.readthedocs.io/en/latest/tutorial.html) to help get you started.

## OpenCAP (command line version)
Compiling OpenCAP requires installing all of our [dependencies](#Dependencies). Assuming these 
are all installed in locations visible to CMake, installation can proceed as follows:

### Basic steps
These steps have been tested on MacOS 12.3, and on Centos7.

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

Once the Makefile is generated, build and install the executable.

```
make

make install

```

### MacOS/Linux 

Any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler. 

For Linux users, depending on your distribution, you may also need to install 
Python3 development libraries e.g. `sudo apt-get install python3.x-dev`.

### Verify installation

To verify your installation, run the tests in the build directory

    make test

To run OpenCAP, specify an input file from the command line

    /path/to/opencap input_file.in 

See the [examples](https://github.com/gayverjr/opencap/tree/main/examples/opencap) directory 
and our [documentation](https://gayverjr.github.io/opencap/input.html) to help get you started.


# Dependencies

Building OpenCAP/PyOpenCAP from source requires working installations of the following:

* C++ compiler with full C++17 language support and standard libraries (**Warning: Default Apple Clang on MacOS is not fully supported**)

* Python3 interpreter and development libraries: version >= 3.6

* [CMake](https://cmake.org/):  version >= 3.12

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/): hierarchical data format, version >= 1.10

* [Eigen](http://eigen.tuxfamily.org/dox/): linear algebra library, version >= 3.3

All of these dependencies are available through standard package managers such as 
[Homebrew](https://brew.sh/), [Conda](https://docs.conda.io/en/latest/), and yum/apt-get on Linux. 

The following packages are automatically built at the CMake step (no action required):

* [Numgrid](https://github.com/dftlibs/numgrid): numerical integration library

* [h5pp](https://github.com/DavidAce/h5pp): C++17 wrapper for HDF5

* [pybind11](https://github.com/pybind/pybind11): C++ Python bindings

The following is included as part of the source distribution (no action required):

* [Date](https://github.com/HowardHinnant/date): Howard Hinnant's date and time library based on the C++11/14/17 <chrono> header

* [asa239](https://people.sc.fsu.edu/~jburkardt/cpp_src/asa239/asa239.html) C++ library which evaluates incomplete Gamma function from John Burkardt. 
  Distributed under GNU [LGPL](https://github.com/gayverjr/opencap/blob/main/opencap/asa/LICENSE) license.
 
# License
OpenCAP and PyOpenCAP are distributed under the [MIT](https://github.com/gayverjr/opencap/blob/main/LICENSE) license. 

[Numgrid](https://github.com/dftlibs/numgrid) is distributed under the Mozilla Public License 2.0. The 
source code is available on [GitHub](https://github.com/dftlibs/numgrid).

# Acknowledgements
This project is funded by the Molecular Sciences Software Institute.

![MolSSI logo](https://github.com/gayverjr/OpenCAP/blob/main/images/molssi_logo.png)

We would like to thank Professor Thomas Sommerfeld for his generous contribution of C++ code for analytic box-CAP integrals. Please 
visit his [repository](https://github.com/tsommerfeld/L2-methods_for_resonances) for many implementations of L^2 methods for resonances!

We would like to thank Mushir Thodika from Temple University for his help in developing the interface with Columbus.

We would like to thank Prof. John Burkardt for his implementation of the Incomplete Gamma function. Many of his wonderful open source codes can be found at his [website](https://people.sc.fsu.edu/~jburkardt/).
