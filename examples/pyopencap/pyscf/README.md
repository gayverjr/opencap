# PyOpenCAP Pyscf examples
This example runs CAP/FCI on the ${}^1\Sigma_g^+$ Feshbach resonance of $H_2$.  The aug-cc-pvqz was augmented with a set of 4s4p diffuse functions placed on two ghost atoms. 

## Scripts

- pyscf_example.py: A CAP/FCI calculation for H2.  Generates zeroth order Hamiltonian (diagonal) and CAP matrix for later processing by cap_trajectory.py. 

- cap_trajectory.py: - cap_trajectory.py: A script which generates eigenvalue trajectories from a properly formatted output file. This file can be easily modified to other systems by changing the guess, the reference energy, and the eta_list defined at the top of the file.

## Outputs
- pyscf_example.out: Output of pyscf_example.py script.

## Data files
- psi4.bas: Basis set formatted in Psi4 style. Used for OpenCAP 
- nwchem.bas: Basis set formatted in NWChem style. Used for Pyscf

## Viewing the trajectory
You can use the supplied script cap_trajectory.py as a template for generating the eigenvalue trajectories and extracting the optimal value of the CAP strength parameter.

    python cap_trajectory.py pyscf_example.out

![trajectory](https://github.com/gayverjr/opencap/blob/master/examples/pyopencap/pyscf/resonance_trajectory.png)
