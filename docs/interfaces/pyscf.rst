PySCF
=======================

PySCF_ is an ab initio computational chemistry program natively implemented in Python. The major
advantage of using Pyscf in tandem with OpenCAP is that calculations can be performed in 
one-shot within the same python script. Since PySCF allows direct control over data structures such as density matrices, the interface between 
PySCF and OpenCAP is seamless. Currently, only FCI has been benchmarked, and here we 
outline how to perform a calculation using this module.

.. _PySCF: http://pyscf.org/

Preliminary: Running the PySCF calculation
------------------------------------------
Please consult the PySCF documentation_ for how run calculations with PySCF. An example 
script_ using FCI is provided in our repository. For FCI, the 
zeroth order Hamiltonian is a diagonal matrix whose entries are the energies of the FCI states.

.. _script: https://github.com/gayverjr/opencap/blob/master/examples/pyopencap/pyscf/pyscf_example.py
.. _documentation: https://sunqm.github.io/pyscf/

Step 1: Defining the System object
----------------------------------

**Molden(recommended)**

The best way to construct the :class:`~pyopencap.System` object is to import the geometry 
and basis set from molden.

.. code-block:: python

    molden_dict = {"basis_file":"molden_in.molden","molecule": "molden"}
    pyscf.tools.molden.from_scf(myhf,"molden_in.molden")
    s = pyopencap.System(molden_dict)

**Inline**

The molecule and basis set can also be specified inline. The "molecule" keyword must 
be set to "read", and then an additional keyword "geometry" must
be specified, with a string that contains the geometry in xyz format. The "basis_file" keyword 
must be set to a path to a basis set file formatted in Psi4 style, which can be downloaded from
the MolSSI BSE_. Other optional keyword for this section include "bohr_coordinates" and
"cart_bf". Please see the :ref:`keywords <keywords>` section for more details. It is recommended to check the
overlap matrix to ensure that the ordering and normalization matches. Up to G-type functions are supported.

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
    
.. _BSE: https://www.basissetexchange.org/

Step 1: Defining the CAP object
-----------------------------------------

The CAP matrix is computed by the :class:`~pyopencap.CAP` object. The constructor 
requires a :class:`~pyopencap.System` object, a dictionary containing the CAP parameters, the number of states,
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
    pc = pycap.CAP(my_system,cap_dict,10,"pyscf")
    
Step 2: Passing the density matrices
------------------------------------
For FCI and related modules, transition densities can be obtained using the :func:`~pyscf.fci.direct_spin1.FCISolver.trans_rdm1`
function of the FCI_ module:

.. code-block:: python
	
	fs = fci.FCI(mol, myhf.mo_coeff)
	e, c = fs.kernel()
	# tdm between ground and 1st excited states
	dm1 = fs.trans_rdm1(fs.ci[0],fs.ci[1],myhf.mo_coeff.shape[1],mol.nelec)
	
.. _FCI: https://sunqm.github.io/pyscf/fci.html

Importantly, trans_rdm1 returns the density matrix in **MO basis**. Thus before passing it to 
PyOpenCAP, it **must be transformed into AO basis**:

.. code-block:: python

    dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
    
Densities are loaded in one at a time using :func:`~pyopencap.CAP.add_tdm`. 
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
using the :func:`~pyopencap.CAP.compute_perturb_cap` function. The matrix can be retrieved using the
:func:`~pyopencap.CAP.get_perturb_cap` function.

.. code-block:: python

    pc.compute_perturb_cap()
    W_mat=pc.get_perturb_cap()
    
*Note:*

When using cartesian d, f, or g-type basis functions, special care must be taken to ensure that the normalization 
conventions match what is used by OpenMolcas. In these cases, :func:`~pyopencap.CAP.compute_ao_cap` 
and then :func:`~pyopencap.CAP.renormalize` or :func:`~pyopencap.CAP.renormalize_cap` 
should be invoked before calling :func:`~pyopencap.CAP.compute_perturb_cap`.

.. code-block:: python

    pc.compute_ao_cap()
    pc.renormalize_cap(pyscf_smat,"pyscf")
    pc.compute_perturb_cap()

Step 4: Generate eigenvalue trajectories
----------------------------------------
Extracting resonance position and width requires analysis of the eigenvalue trajectories. 
A template trajectory analysis script is provided in the 
repository_. Development of automated tools for trajectory analysis is a subject of future work.

.. _repository: https://github.com/gayverjr/opencap/blob/master/examples/pyopencap/pyscf/cap_trajectory.py


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