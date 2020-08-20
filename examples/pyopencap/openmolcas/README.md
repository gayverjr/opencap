# PyOpenCAP Molcas examples
The following examples run CAP/(X)MS-CASPT2 calculations on the pi* shape resonance of N2-. The basis set is aug-cc-pvtz + a set of even tempered 3s3p3d diffuse basis functions placed on a ghost atom in the center of the molecule.

## Scripts

- example.py: A CAP/XMS-CASPT2(5,11) calculation with symmetry turned off. Plots uncorrected and corrected trajectories, and computes their stationary points.

- symmetry_example.py: A CAP/MS-CASPT2(5,11) calculation using D2h symmetry. Plots uncorrected and corrected trajectories, and computes their stationary points.

- cartesian_example.py: A CAP/MS-CASPT2(5,11) calculation using cartesian basis functions. Special care is required in this case to properly re-normalize the CAP matrix to be consistent with the electronic structure package. Plots uncorrected and corrected trajectories, and computes their stationary points.

## Notebook

- Tutorial.ipynb: A tutorial to help get you started on using OpenCAP with OpenMolcas.

