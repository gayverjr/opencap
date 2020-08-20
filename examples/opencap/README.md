# OpenCAP examples

The following examples are for running CAP/(X)MS-CASPT2 calculations on the pi* shape resonance of N2-. The basis set is aug-cc-pvtz + a set of even tempered 3s3p3d diffuse basis functions placed on a ghost atom in the center of the molecule.

## Input files

- symmetry_on.in: CAP/MS-CASPT2(5,11) job with symmetry, Box CAP

- symmetry_off.in: CAP/XMS-CASPT2(5,11) job with symmetry turned off, Box CAP

- cartesian.in: CAP/MS-CASPT2(5,11) job with Cartesian basis functions, Voronoi CAP

## Output files

- ms-caspt2.out: OpenCAP output for symmetry_on.in input file

- xms-caspt2.out: OpenCAP output for symmetry_off.in input file.

- cartesian_example.out: OpenCAP output for cartesian.in input file.

## Data files

- cart.out: MS-CASPT2(5,11) molcas output file using cartesian basis functions

- nosymm.out: XMS-CASPT2(5,11) molcas output file, symmetry turned off

- nosymm.rasscf.molden: Molden input file 

- symm.out: MS-CASPT2(5,11) molcas output file using D2h symmetry

- symm.rassi.h5: RASSI data file for MS-CASPT2(5,11) using D2h symmetry

- cart.rassi.h5: RASSI data file for MS-CASPT2(5,11) using cartesian basis functions

- nosymm.rassi.h5: RASSI data file for XMS-CASPT2(5,11), symmetry turned off

- basis.bas: Basis set file in Psi4 format

## Scripts
- cap_trajectory.py: A script which generates eigenvalue trajectories from an OpenCAP output file. This file can be easily modified to other systems by changing the guess, the reference energy, and the eta_list defined at the top of the file.

## Viewing the trajectories
You can use the supplied script cap_trajectory.py as a template for generating the eigenvalue trajectories and extracting the optimal value of the CAP strength parameter.

    python cap_trajectory.py xms-caspt2.out

![trajectory](https://github.com/gayverjr/opencap/blob/master/examples/pyopencap/openmolcas/Images/res_trajectory.png)
