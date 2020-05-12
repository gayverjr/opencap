# OpenMolcas Interface

OpenMolcas is a widely used quantum chemistry package which specializes in multiconfigurational approaches. Here, we focus on how to use the software in tandem with OpenCAP to perform calculations on resonances. 

The key to the OpenMolcas interface is the RAS State Interaction program (RASSI). As stated in the OpenMolcas reference, RASSI "computes matrix elements of the Hamiltonian and other operators in a wave function basis, which consists of individually optimized CI expansions from the RASSCF program". Of crucial importance to us is that RASSI is able to compute and export one particle transition density matrices between each pair of states. When the 'TRD1' keyword is activated, these densities are exported to an HDF5 file titled $JOBNAME.rassi.h5. OpenCAP uses these densities to compute the %CAP matrix elements in wave-function basis, which is required to perform projected %CAP calculations.

In addition to the one particle densities, OpenCAP also requires a zeroth order Hamiltonian. For a proper description of correlation, an effective Hamiltonian from a multi-state second order perturbation theory method is required. OpenCAP is capable of reading the effective Hamiltonian from OpenMolcas outputs for the following methods: MS-CASPT2,XMS-CASPT2.

There is a working example in the examples directory of the repository, which includes OpenMolcas and OpenCAP inputs and outputs. 

#### OpenMolcas example files ####

| Keyword | Description |
|----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| anion_reference.out | Output from 10 state MS-CASPT2 OpenMolcas calculation on N2-.  Contains MS-CASPT2 effective Hamiltonian which can be read in by OpenCAP. |
| test.rassi.h5 | RASSI HDF5 file containing TDMs.  |
| ground_state_reference.out | Output from 10 state MS-CASPT2 OpenMolcas calculation on neutral N2. Useful for obtaining ground state energy, which is required to compute resonance position. |
| test.in | Test OpenCAP input which utilizes OpenMolcas interface. |
| test_bas.bas | Basis set file in Psi4 format. Aug-cc-pvtz +3s3p3d set of diffuse functions on ghost atom. |
