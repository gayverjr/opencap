.. _CAPHamiltonian:

CAPHamiltonian
==============================
This section briefly describes how to use the CAPHamiltonian object to generate 
eigenvalue trajectories.

Initialization
----------------
:class:`~pyopencap.analysis.CAPHamiltonian` objects can be initialized in one of two ways. 
The first is to pass H0 and W as numpy arrays:

.. code-block:: python
	
	from pyopencap.analysis.CAPTrajectory import CAPHamiltonian	
	CAPH = CAPHamiltonian(H0=h0,W=mat)

The other is to read them in from an OpenCAP output file, or from a Q-Chem output file generated 
by a Projected CAP-EOM-CC calculation. 

.. code-block:: python

	CAPH = CAPHamiltonian(output_file="path/to/output.out")

If one wishes to exclude some of the states from the analysis, this can be accomplished through the 
by placing their indices in a list (starting from 0) and passing it into the 
`exclude_states` keyword argument:

.. code-block:: python

	exclude_states = [2,5,7]
	CAPH = CAPHamiltonian(H0=h0,W=mat,exclude_states=exclude_states)

Importantly, in all cases, **the W matrix is assumed to be pre-multiplied by a factor of -1.0**.

Diagonalization 
----------------
The :func:`~pyopencap.analysis.CAPHamiltonian.run_trajectory` function diagonalizes the 
CAP Hamiltonian over a range of eta values (and at a specified value of the cap lambda parameter
if using a CR-CAP).

.. code-block:: python
	
	eta_list = np.linspace(0,2000,101)
	eta_list = eta_list * 1E-5
	CAPH.run_trajectory(eta_list,cap_lambda=0.0)

Since the W matrix is assumed to be multiplied by a factor of -1.0 
upon instantiation, the following matrix is actually diagonalized at each step: 

.. centered:: :math:`H^{CAP}=H_0+(i\eta- \lambda) W`

and each eigenpair is stored in a :class:`~pyopencap.analysis.Root` object.
After all of the diagonalizations are finished, individual states can be tracked 
using the :func:`~pyopencap.analysis.CAPHamiltonian.track_state` function:

.. code-block:: python

	traj = CAPH.track_state(4,tracking="overlap")

The `traj` variable is a :class:`~pyopencap.analysis.EigenvalueTrajectory` object, which 
contains helpful functions for analysis. Indices for states start from 0, and there are 
two options for tracking states: "overlap" (the default), and "energy". See :class:`~pyopencap.analysis.EigenvalueTrajectory`
for more details.

Visualization
--------------
The energies of all states computed are stored in the `all_energies` class attribute. This can 
very useful for graphical searches e.g.

.. code-block:: python
	
	import matplotlib.pyplot as plt
	import numpy as np
	plt.plot(np.real(CAPH.total_energies),np.imag(CAPH.total_energies),'ro')
	plt.show()

There is also a function :func:`~pyopencap.analysis.CAPHamiltonian.energies_ev` which returns 
the excitation energies in eV with respect to specified reference energy.

.. code-block:: python

	E_ev = CAPH.energies_ev(ref_energy)
	plt.plot(np.real(E_ev),np.imag(CAPH.energies_ev(E_ev),'ro')
	plt.show()

.. autoclass:: pyopencap.analysis.CAPHamiltonian
   :members:

   .. automethod:: __init__
