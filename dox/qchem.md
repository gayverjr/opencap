Q-Chem Interface {#qchem}
=============================

[Q-Chem](https://molcas.gitlab.io/OpenMolcas/sphinx/) is a widely used quantum chemistry package with a variety of open-shell and excited state electronic structure methods. We officially support an interface with EOM-CC calculations performed by Q-Chem. There is also an experimental TDDFT interface. 

System
======

OpenCAP can extract the basis set from .fchk files generated by Q-Chem.

__Relevant keywords__
| Keyword    | Valid options                     | Description                                                                                                                                                                                                                                                                                      |
|------------|-----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| molecule   | molden,rassi_h5,inline,qchem_fchk | Specifies which format to read the molecular geometry. If "inline" is chosen, the "$geometry" section is also required.                                                                                                                                                                          |
| basis_file | path to basis file                | Specifies the path to the basis file. When "molecule" is set to "molden","rassi_h5", or "qchem_fchk", this keyword should be set to a path to a file of the specified type. When "molecule" is set to  "inline", this keyword should be set to a path to a basis set file formatted in "Psi4" style. |

Projected CAP
=============

The keys to the Q-Chem interface are the GUI, CC_TRANS_PROP, and STATE_ANALYSIS keywords.  

     GUI  2
	 CC_TRANS_PROP 2
	 STATE_ANALYSIS 1

Including these keywords in the $rem card will export all 1RDMs and 1TDMs between all EOM states to a .fchk file. 

Densities can be formatted in one of two ways (automatically detected by OpenCAP):

1. Triangular (Format of CC_TRANS_PROP=2, %CAP matrix is assumed symmetric in this case)

    State Density 1...
	
    State Density N
	
	TDM 1-->2
	
	TDM 1-->3...
	
	TDM 1-->N...
	
	TDM 2-->3...
	
	TDM 2-->N...

2. All densities (experimental format)

    State Density 1 ...

    State Density N 

	TDM 1-->2
	
	TDM 1-->3...
	
	TDM 1-->N
	
	TDM 2-->1
	
	TDM 2-->3...
	
	TDM 2-->N...



Separate alpha/beta densities and spin traced one particle densities are both supported, and 
are automatically detected by OpenCAP.

__Relevant keywords__
| Keyword | Type | Description |
|---------------|--------|---------------------------------------------------------|
| qchem_fchk | string | Relative or absolute path to .fchk file containing TDMs.     |
| qchem_output | string | Relative or absolute path to Q-Chem output file.   |
| h0_file | string | Relative or absolute path to properly formatted zeroth order Hamiltonian file.   |


__Example input__

    $system
    molecule qchem_fchk
    basis_file qc_input.fchk
    $end

    $projected_cap
    method eom
    package qchem
    qchem_fchk qc_input.fchk
    nstates 10
    CAP_TYPE voronoi
    R_CUT 3.00
    Radial_precision  14
    Angular_points    110
    $end
    