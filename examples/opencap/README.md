# OpenCAP examples
All of these examples are for N2-.

## Input files
- symmetry_on.in: MS-CASPT2 job with symmetry, Box CAP
- symmetry_off.in: XMS-CASPT2 job with symmetry turned off, Box CAP
- cartesian.in: MS-CASPT2 job with Cartesian basis functions, Voronoi CAP

## Output files
- ms-caspt2.out: OpenCAP output for symmetry_on.in input file

## Data files
- cart.out: MS-CASPT2 molcas output file
- nosymm.out: XMS-CASPT2 molcas output file
- nosymm.rasscf.molden: Molden input file for XMS-CASPT2
- symm.out: MS-CASPT2 molcas output file
- symm.rassi.h5: RASSI data file for MS-CASPT2
- cart.rassi.h5: RASSI data file for MS-CASPT2
- nosymm.rassi.h5: RASSI data file for XMS-CASPT2
- basis.bas: Basis set file in Psi4 format

## Viewing the trajectory
Use the supplied script cap_trajectory.py as template for generating the eigenvalue trajectories and extracting the optimal value of the CAP strength parameter.

    cap_trajectory.py ms-caspt2.out