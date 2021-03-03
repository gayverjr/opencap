.. _analysis:

Analysis Tools
=================================
PyOpenCAP provides user friendly tools for analysis of eigenvalue trajectories. This page
is currently under construction, but the basics are as follows:

**Basic usage**

The :class:`~pyopencap.analysis.CAPHamiltonian` contains functions aimed at diagonalization 
of the CAP Hamiltonian over a range of eta values. Assuming one has already obtained H0 and 
W in the state basis as numpy matrices:

.. code-block:: python
	
	from pyopencap.analysis.CAPTrajectory import CAPHamiltonian	
	eta_list = np.linspace(0,2000,101)
	eta_list = eta_list * 1E-5
	CAPH = CAPHamiltonian(H0=h0,W=mat)
	CAPH.run_trajectory(eta_list,cap_lambda=0.0)
	# track the 4th state
	traj = CAPH.track_state(4,tracking="overlap")

Alternatively, one can read in H0 and W from OpenCAP/Q-Chem output files:

.. code-block:: python

	from pyopencap.analysis.CAPTrajectory import CAPHamiltonian	
	eta_list = np.linspace(0,2000,101)
	eta_list = eta_list * 1E-5
	CAPH = CAPHamiltonian(output_file="path/to/output.out")
	CAPH.run_trajectory(eta_list,cap_lambda=0.0)
	# track the 4th state
	traj = CAPH.track_state(4,tracking="overlap")
	
In both snippets, `traj` is now a :class:`~pyopencap.analysis.EigenvalueTrajectory` object, which 
contains helpful functions for analysis. For example, one can find the optimal value of the CAP 
strength parameter for uncorrected/corrected trajectories:

.. code-block:: python

	uc_energy,uc_eta_opt = traj.find_eta_opt()
	corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

For more information, please see the documentation for the :class:`~pyopencap.analysis.CAPHamiltonian`
and :class:`~pyopencap.analysis.EigenvalueTrajectory` classes.

.. toctree::
   :maxdepth: 1

   CAPHamiltonian
   EigenvalueTrajectory
   Root