Installation
=============

Install with pip (recommended)
---------------------------------

.. code:: bash
   
    pip install pyopencap

Precompiled Python wheels are available on Pypi for almost all Linux systems and 
most MacOS systems, for Python versions 3.5 and later.

Build from source
-----------------

**Dependencies**

Compiling PyOpenCAP from source requires first installing the following dependencies:

* C++ compiler with full C++17 language support and standard libraries (**Warning: Default Apple Clang on MacOS is not supported**)

* Python3 version >= 3.4

* CMake_: version >= 3.12

* HDF5_: hierarchical data format, version >= 1.10

* Eigen_: linear algebra library, version >= 3.3

All of these dependencies are available through standard package managers such as 
Homebrew_, Conda_, and yum/apt-get 
on Linux. Please consult their respective websites for more information, and ensure that you 
have recent enough versions installed on your machine. 

.. _Homebrew: https://brew.sh/
.. _Conda: https://docs.conda.io/en/latest/ 
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _CMake: https://cmake.org/
.. _Eigen: http://eigen.tuxfamily.org/dox/

**Compiler**

For Linux users, any compiler which fully supports the C++17 standard should work 
(e.g GCC 7.x or later). If you are unsure, try updating to the latest version of your 
compiler.

For Mac users, as of MacOS 10.15 Catalina, the Apple Clang provided by XCode will not work due to missing standard 
library features. We suggest installing the latest version of GCC (currently 10.2) 
from Homebrew_, and then setting the following 
environment variables before attempting to build from source:

.. code:: bash

   # for GCC 10 installed by brew

   export CC=gcc-10

   export CXX=g++-10


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
and running :code:`pytest`. The following python packages are required to run the tests:

.. code:: bash

   pip install h5py
   pip install numpy
   pip install pytest
   pip install pyscf