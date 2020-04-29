
# Welcome to OpenCAP!

OpenCAP is an open-source application aimed at extending the capabilities of electronic structure packages to describe metastable electronic states. We currently support an interface with the [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas) package to compute resonance positions and widths using the complex absorbing potential method (CAP). See the example below for a brief tutorial on how to use the software.





# Installation

## Dependencies

OpenCAP requires the following:

*  C++ compiler with C++17 support

* [CMake](https://cmake.org/)  version 3.10 or higher

*  [HDF5](https://www.hdfgroup.org/solutions/hdf5/) hierarchical data format

* [Armadillo](http://arma.sourceforge.net/) linear algebra library configured with HDF5 support

* [Numgrid](https://github.com/dftlibs/numgrid) numerical integration library



## Install with CMake

First clone the git repo

```

git clone https://github.com/gayverjr/OpenCAP.git

```

Generate the Makefile using CMake

```

cd OpenCAP

mkdir build

cd build

cmake ../

```

Finally, make and install

```

make

make install

```



## Troubleshooting

If Armadillo is not installed in a standard location visible by CMake (e.g. on a computing cluster or machine without root access), the ARMADILLO_INCLUDE_DIR and ARMADILLO_LIBRARY variables must be passed to CMake.

```

cmake -DARMADILLO_INCLUDE_DIR=/path/to/arma/include -DARMADILLO_LIBRARY=/path/to/arma/lib/libarmadillo.so \

../

```

Additionally, Armadillo must be configured with HDF5 support, which may not be the case when installed by package managers such as Anaconda. Therefore, we suggest compiling Armadillo from source and configuring with HDF5 support. See their [documentation](http://arma.sourceforge.net/docs.html) for how to do so.

# Example

Below, we outline how to conduct a CAP-augmented multi-state complete active space second order perturbation theory (CAP-MS-CASPT2) calculation on the <sup>2</sup>&#928<sub>g</sub> resonance state of N<sub>2</sub><sup>-</sup> using OpenCAP in tandem with OpenMolcas.



* Step 1:  Perform an MS-CASPT2 calculation on the anion in OpenMolcas, using the ``TRD1`` keyword in the ``RASSI`` module (see [example/anion_reference.out](https://github.com/gayverjr/OpenCAP/blob/master/example/anion_reference.out)), which exports state one particle density matrices and transition density matrices to a file called $JOBNAME.rassi.h5.

*  Step 2:  Construct an input file for OpenCAP(see [example/test.in](https://github.com/gayverjr/OpenCAP/blob/master/example/test.in)) which contains the relevant system information, and points to the rassi.h5 and OpenMolcas output file which contains the  MS-CASPT2 effective Hamiltonian.

* Step 3: Run OpenCAP to generate the matrices required for the projected CAP calculation.

```OpenCAP example/test.in > example/test.out ```

* Step 4: Use the script provided in [example/cap_trajectory.py](https://github.com/gayverjr/OpenCAP/blob/master/example/cap_trajectory.py) to analyze the eigenvalue trajectory to extract resonance position and width.

```` python example/cap_trajectory.py example/test.out````

The script provides visualization and numerical analysis functionality to identify the resonance trajectory and extract the resonance parameters. It should be modified to suit your purposes.  For our MS-CASPT2 example, the corrected (blue) and uncorrected trajectories look like:

![CAP trajectory](https://github.com/gayverjr/OpenCAP/blob/master/images/cap_traj.png)

Analysis of the uncorrected trajectory results in a resonance energy of 2.26eV, and a width of 0.338eV, which is in good agreement with previous theoretical results obtained for this system.

# Acknowledgements
This project was funded by the Molecular Sciences Software Institute.


![MolSSI logo](https://github.com/gayverjr/OpenCAP/blob/master/images/molssi_logo.png)
