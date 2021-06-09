.. _getting_started:

Getting Started
=======================================


**Constructing the System object** 

The :class:`~pyopencap.System` object of PyOpenCAP contains the geometry and basis set information, as well
as the overlap matrix. The constructor takes in a Python dictionary as an argument, 
and understands a specific set of keywords_ . There are four equivalent ways of specifying
the geometry and basis set: qchem_fchk, rassi_h5, molden, and inline. Here, we'll use the rassi_h5 file.

.. code-block:: python

    sys_dict = {"molecule": "rassi_h5","basis_file": "path/to/rassi/file.h5"}
    s = pyopencap.System(sys_dict)
    smat = s.get_overlap_mat()
    
.. _keywords: https://gayverjropencap.readthedocs.io/en/latest/keywords.html

**Constructing the CAP object**

The CAP matrix is computed by the :class:`~pyopencap.CAP` object. The constructor 
requires a :class:`~pyopencap.System` object, a dictionary containing the CAP parameters, 
and the number of states.

.. code-block:: python
	
	nstates = 10
    cap_dict = {"cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"}
    pc = pyopencap.CAP(s,cap_dict,nstates)
    
**Parsing electronic structure data from file**

The :func:`~pyopencap.CAP.read_data` function can read in the effective Hamiltonian
and densities in one-shot when passed a Python dictionary with the right keywords_. 

.. code-block:: python

    es_dict = {"method" : "xms-caspt2",
              "package": "openmolcas",
      "molcas_output":"path/to/output/file.out",
      "rassi_h5":  "path/to/rassi/file.h5"}
    pc.read_data(es_dict)
    h0 = pc.get_H()

**Passing densities in RAM**

Alternatively, one can load in the densities one at a time using the :func:`~pyopencap.CAP.add_tdms` 
or :func:`~pyopencap.CAP.add_tdm` functions. The density matrices should be in atomic orbital basis,
with the same atomic orbital ordering as the :class:`~pyopencap.System` (which can be verify using 
:class:`~pyopencap.check_overlap_matrix`). The example below shows how one might pass the densities 
from a PySCF calculation:

.. code-block:: python
    
    s.check_overlap_mat(pyscf_smat,"pyscf")
    pc = pyopencap.CAP(s,cap_dict,10)
    for i in range(0,10):
        for j in range(i,10):
            dm1 = fs.trans_rdm1(fs.ci[i],fs.ci[j],myhf.mo_coeff.shape[1],mol.nelec)
            dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
            pc.add_tdm(dm1_ao,i,j,"pyscf")
            if i!=j:
                pc.add_tdm(dm1_ao,j,i,"pyscf")
    

Once all of the densities are loaded, the CAP matrix is computed 
using the :func:`~pyopencap.CAP.compute_projected_cap` function. The matrix can be retrieved using the
:func:`~pyopencap.CAP.get_projected_cap` function.

    >>> pc.compute_projected_cap()
    >>> W_mat=pc.get_projected_cap()

We now have our zeroth order Hamiltonian (stored in h0) and our CAP matrix(W_mat) in
the state basis. Extracting resonance position and width requires analysis of the 
eigenvalue trajectories. 

**Analysis**

PyOpenCAP provides user friendly tools for analysis of eigenvalue trajectories.

The :class:`~pyopencap.analysis.CAPHamiltonian` contains functions aimed at diagonalization 
of the CAP Hamiltonian over a range of eta values. Assuming one has already obtained H0 and 
W in the state basis as numpy matrices, it can be constructed as such:

.. code-block:: python
	
	from pyopencap.analysis.CAPTrajectory import CAPHamiltonian	
	eta_list = np.linspace(0,2000,101)
	eta_list = eta_list * 1E-5
	CAPH = CAPHamiltonian(H0=h0,W=mat)
	# equivalently
	CAPH = CAPHamiltonian(pc=pc)
	CAPH.run_trajectory(eta_list,cap_lambda=0.0)
	# track the 4th state
	traj = CAPH.track_state(4,tracking="overlap")
	
`traj` is now a :class:`~pyopencap.analysis.EigenvalueTrajectory` object, which 
contains helpful functions for analysis. For example, one can find the optimal value of the CAP 
strength parameter (and therefore, best estimate of resonance position and width) for uncorrected/corrected trajectories:

.. code-block:: python

	uc_energy,uc_eta_opt = traj.find_eta_opt()
	corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

For more information, please see the documentation for the :class:`~pyopencap.analysis.CAPHamiltonian`
and :class:`~pyopencap.analysis.EigenvalueTrajectory` classes.

**See more**

Please see the notebooks in our repository_ for detailed examples which demonstrate the full 
functionality of PyOpenCAP.

.. _repository: https://github.com/gayverjr/opencap/tree/main/examples/analysis




