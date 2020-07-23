PySCF
=======================

PySCF is an ab initio computational chemistry program natively implemented in Python. The major
advantage of using Pyscf in tandem with OpenCAP is that calculations can be performed in 
one-shot within the same python script, without the need to perform any I/O. Since PySCF 
allows direct control over data structures such as density matrices, the interface between 
PySCF and OpenCAP is seamless, and there are numerous possibilities. Currently, only FCI
has been benchmarked, and here we outline how to perform a calculation using this module.

Preliminary: Running the PySCF calculation
------------------------------------------
Please consult the PySCF documentation for how run calculations with PySCF. An example 
script using FCI is provided in our repository at examples/pyscf/h2.py. For FCI, the 
zeroth order Hamiltonian is a diagonal matrix whose entries are the energies of the FCI states.


Step 1: Defining the System and Projected_CAP objects
-----------------------------------------------------

**System object**

The molecule and basis set must be specified manually. The "molecule" keyword must 
be set to "read", and then an additional keyword "geometry:" must
be specified, with a string that contains the geometry in xyz format. The "basis_file" keyword 
must be set to a path to a basis set file formatted in Psi4 style, which can be downloaded from
the basis set exchange. Other optional keyword for this section include "bohr_coordinates" and
cart_bf. Please see the keywords section for more details. It is recommended to check the
overlap matrix to ensure that the ordering matches.

.. code-block:: python
	
    pyscf_smat = scf.hf.get_ovlp(mol)
    sys_dict = {"geometry":    '''N  0  0   1.039
                              N  0  0   -1.039
                              X   0  0   0.0''',
            		"molecule" : "read",
            		"basis_file":"path/to/basis.bas",
            		"cart_bf":"d",
            		"bohr_coordinates:": "true"}
    s.check_overlap_mat(pyscf_smat,"pyscf")
    
**Projected_CAP object**

The CAP matrix is computed by the "Projected_CAP" object. The constructor for the Projected_CAP 
object requires a System object, a dictionary containing the CAP parameters, the number of states,
and finally the string "pyscf", which denotes the ordering of the atomic orbital basis
set. An example is provided below. Please see the keywords section for more information on
the CAP parameters.

.. code-block:: python

    cap_dict = {"cap_type": "box",
            	"cap_x":"2.76",
            	"cap_y":"2.76",
            	"cap_z":"4.88",
            	"Radial_precision": "14",
            	"angular_points": "110"}
    pc = pycap.Projected_CAP(my_system,cap_dict,10,"pyscf")
    
Step 2: Passing the density matrices
------------------------------------
For FCI and related modules, transition densities can be obtained using the trans_rdm function:

.. code-block:: python
	
	fs = fci.FCI(mol, myhf.mo_coeff)
	e, c = fs.kernel()
	# tdm between ground and 1st excited states
	dm1 = fs.trans_rdm1(fs.ci[0],fs.ci[1],myhf.mo_coeff.shape[1],mol.nelec)

Importantly, trans_rdm1 returns the density matrix in **MO basis**. Thus before passing it to 
PyOpenCAP, it **must be transformed into AO basis**:

.. code-block:: python

    dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
    
Densities are loaded into the Projected_CAP object one at a time using the "add_tdm" function. 
Ensure that the indices of each state match those of the zeroth order Hamiltonian.

.. code-block:: python

    for i in range(0,len(fs.ci)):
        for j in range(0,len(fs.ci)):
            dm1 = fs.trans_rdm1(fs.ci[i],fs.ci[j],myhf.mo_coeff.shape[1],mol.nelec)
            dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
            pc.add_tdm(dm1_ao,i,j,"pyscf")


Step 3: Computing the CAP matrix
--------------------------------
Once all of the densities are loaded, the CAP matrix is computed 
using the "compute_projected_cap" function. The matrix can be retrieved using the
"get_projected_cap" function.

.. code-block:: python

    pc.compute_projected_cap()
    W_mat=pc.get_projected_cap()

Step 4: Generate eigenvalue trajectories
----------------------------------------
Eigenvalue trajectories by diagonalizing H0-inW over a range of eta values, and then 
analyzing the results. A template (pyscf_ex.py) script is provided in the 
examples/pyscf directory of our repository. Automated tools for trajectory analysis
is a subject of future work.


Officially supported methods
----------------------------

* Full CI

Coming (hopefully) soon
-----------------------

* EOM-CCSD
* ADC

Untested (use at your own risk!)
--------------------------------
Any module which one particle transition densities available can be supported. 
This includes all methods which can utilize the trans_rdm1 function, including but not limited to:

* MRPT
* CISD
* TDDFT