# PyOpenCAP Pyscf examples

## Scripts
- pyscf_example.py: A CAP/FCI calculation for H2. 
- cap_trajectory.py: A script for analyzing the eigenvalue trajectories.

## Outputs
- pyscf_example.out

## Data files
- psi4.bas: Basis set formatted in Psi4 style. Used for OpenCAP 
- nwchem.bas: Basis set formatted in NWChem style. Used for Pyscf

## Viewing the trajectory
You can use the supplied script cap_trajectory.py as a template for generating the eigenvalue trajectories and extracting the optimal value of the CAP strength parameter.

    python cap_trajectory.py pyscf_example.out