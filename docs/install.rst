Installation
=============

Install with pip (recommended)
---------------------------------

.. code:: bash
   
    pip install pyopencap
    # or
    pip3 install pyopencap

Precompiled Python wheels are available on Pypi for almost all Linux systems and 
most MacOS systems, for Python versions 3.6 and later.

Build from source
-----------------

**Dependencies**

Compiling PyOpenCAP from source requires first installing the following dependencies:

* C++ compiler with full C++17 language support and standard libraries

* Python3 interpreter and development libraries: version >= 3.6

* CMake_: version >= 3.12

* HDF5_: hierarchical data format, version >= 1.10

* Eigen_: linear algebra library, version >= 3.3

All of these dependencies are available through standard package managers such as 
Homebrew_, Conda_, and yum/apt-get on Linux. 

.. _Homebrew: https://brew.sh/
.. _Conda: https://docs.conda.io/en/latest/ 
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _CMake: https://cmake.org/
.. _Eigen: http://eigen.tuxfamily.org/dox/

**Compiler**

For Mac/Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler. 

**Building the package**

If your operating system/Python environment is not covered by any of our pre-built wheels,
the command :code:`pip install pyopencap` will download the tarball from Pypi and try to compile from source.
You can also clone the repository and install a local version:

.. code:: bash

    git clone https://github.com/gayverjr/opencap.git

    cd opencap

    pip install .

Compiling from source will take several minutes. To monitor your progress, you can run pip 
with the `--verbose` flag.

To ensure that the installation was successful, return to your home directory, start a Python shell, and type:

.. code:: python
   
   import pyopencap
    
If you cloned the repository, you can run the tests by entering the `pyopencap` directory, 
and running :code:`pytest`.