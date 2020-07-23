<div align="left">
  <img src="https://github.com/gayverjr/opencap/blob/fchk_parse/images/opencap_logo.png" height="200px">
</div>

# OpenCAP: An open-source program for studying resonances in molecules.                                                               
| **Documentation** | [![][docs-img]][docs-url] [![][rtd-img]][rtd-url] |
| :------ | :------- |
| **Build Status**  | [![][travis-img]][travis-url] [![][lgtm-img]][lgtm-url] |

[docs-img]: https://img.shields.io/badge/doc-latest-blue.svg
[docs-url]: https://gayverjr.github.io/opencap/
[travis-img]: https://travis-ci.com/gayverjr/opencap.svg?branch=master
[travis-url]: https://travis-ci.com/gayverjr/opencap
[lgtm-img]: https://img.shields.io/lgtm/grade/cpp/g/gayverjr/opencap.svg?logo=lgtm&logoWidth=18
[lgtm-url]: https://lgtm.com/projects/g/gayverjr/opencap/context:cpp
[rtd-img]: https://readthedocs.org/projects/gayverjropencap/badge/?version=latest
[rtd-url]: https://gayverjropencap.readthedocs.io/en/latest/?badge=latest

OpenCAP is an open-source application aimed at extending the capabilities of electronic structure packages to describe metastable electronic states. 
We currently support an interface with the [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas) package to 
compute resonance positions and widths using the complex absorbing potential method (CAP). 
See the example below for a brief tutorial on how to use the software.


Documentation (in progress): https://gayverjr.github.io/opencap/





# Installation

## Dependencies

OpenCAP requires the following:

*  C++ compiler with C++17 support

* [CMake](https://cmake.org/)  verison >= 3.12

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/) hierarchical data format, version >= 1.8

* [Eigen](http://eigen.tuxfamily.org/dox/) linear algebra library

* [Numgrid](https://github.com/dftlibs/numgrid) numerical integration library



## Install the executable with CMake

First clone the git repo

```

git clone https://github.com/gayverjr/opencap.git

cd opencap

```

Generate the Makefile using CMake

```

cd opencap

mkdir build

cd build

cmake ../

```

Finally, make and install

```

make

make install

```

## Installing the python package (experimental)

First clone the git repo

```

git clone https://github.com/gayverjr/opencap.git

cd opencap

```

Then pip install

```

pip install .

```

The module is now importable within a python shell as "pyopencap".

```
import pyopencap

```

# Acknowledgements
This project is funded by the Molecular Sciences Software Institute.


![MolSSI logo](https://github.com/gayverjr/OpenCAP/blob/master/images/molssi_logo.png)
