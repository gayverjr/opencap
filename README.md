
# Welcome to OpenCAP!

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

pip install ./opencap

```

The module is now importable within a python shell as "pycap".

```
import pycap

```

# Example

Below, we outline how to conduct a CAP-augmented multi-state complete active space second order perturbation theory (CAP-MS-CASPT2) calculation on the <sup>2</sup>&#928;<sub>g</sub> resonance state of N<sub>2</sub><sup>-</sup> using OpenCAP in tandem with OpenMolcas.



* Step 1:  Perform an MS-CASPT2 calculation on the anion in OpenMolcas, using the ``TRD1`` keyword in the ``RASSI`` module (see [example/anion_reference.out](https://github.com/gayverjr/OpenCAP/blob/master/example/molcas/anion_ms.out)), which exports state one particle density matrices and transition density matrices to a file called $JOBNAME.rassi.h5.

*  Step 2:  Construct an input file for OpenCAP(see [example/test.in](https://github.com/gayverjr/OpenCAP/blob/master/example/molcas/test.in)) which contains the relevant system information, and points to the rassi.h5 and OpenMolcas output file which contains the  MS-CASPT2 effective Hamiltonian.

* Step 3: Run OpenCAP to generate the matrices required for the projected CAP calculation.

```opencap test.in > test.out ```

* Step 4: Use the script provided in [example/cap_trajectory.py](https://github.com/gayverjr/OpenCAP/blob/master/example/molcas/cap_trajectory.py) to analyze the eigenvalue trajectory to extract resonance position and width.

```` python cap_trajectory.py test.out````

The script provides visualization and numerical analysis functionality to identify the resonance trajectory and extract the resonance parameters. It should be modified to suit your purposes.  For our MS-CASPT2 example, the corrected (blue) and uncorrected trajectories look like:

![CAP trajectory](https://github.com/gayverjr/OpenCAP/blob/master/images/cap_traj.png)

Analysis of the uncorrected trajectory results in a resonance energy of 2.16eV, and a width of 0.29eV, which is in good agreement with previous theoretical results obtained for this system.

# Acknowledgements
This project is funded by the Molecular Sciences Software Institute.


![MolSSI logo](https://github.com/gayverjr/OpenCAP/blob/master/images/molssi_logo.png)
