.. _EigenvalueTrajectory:

EigenvalueTrajectory
==============================
This section briefly describes how to use the EigenvalueTrajectory object to analyze 
eigenvalue trajectories.

Initialization
--------------
:class:`~pyopencap.analysis.EigenvalueTrajectory` objects are generated by the :func:`~pyopencap.analysis.CAPHamiltonian.track_state`
function of the :class:`~pyopencap.analysis.CAPHamiltonian` class. The state index i starts from 0, and the 
first state in the trajectory is the ith eigen pair generated by the first diagonalization at :math:`\eta=0`. 

State Tracking
---------------
At each diagonalization of the CAP Hamiltonian, a set of n :math:`\eta`-dependent eigenpairs
is generated, and the next point along the eigenvalue trajectory is chosen from that set 
using one of two criterion: overlap and energy.

When overlap tracking is used (the default), at each step, the state with maximum overlap with the 
previous state is chosen as the next point on the trajectory. The overlap is simply 
calculated as the absolute value of the dot product between the eigenvectors. 

When energy tracking is used, at each step, the state with the minimum difference in energy 
from the previous state is chosen. Since the energies are complex, the difference is computed as 
the modulus of the difference of the two complex energies.

Each state in the trajectory is stored as a :class:`~pyopencap.analysis.Root` object stored in the 
`states` attribute. Lists of uncorrected energies and :math:`\eta` values (in order of 
smallest to largest :math:`\eta` value) are also stored in the `uncorrected_energies` 
and `etas` class attributes for convenience.


Corrected Trajectories
-----------------------
Raw uncorrected energies obtained from diagonalizaiton of the CAP Hamiltonian can be sensitive to 
CAP onset and basis set quality. Practitioners of CAP theory often report so called `corrected` 
energies, the exact form of which may vary from publication to publication. The form of correction 
we utilize here is that of [Cederbaum2002]_, which has the form:

:math:`U(\eta)=E(\eta)-\eta\frac{\partial E(\eta) }{\partial \eta}`.

Corrected trajectories are automatically computed when one obtains an :class:`~pyopencap.analysis.EigenvalueTrajectory`
object generated by the :func:`~pyopencap.analysis.CAPHamiltonian.track_state` funciton, and 
are stored in the class attribute `corrected_energies`.

:math:`\eta_{opt}`
------------------
As briefly outlined in the :ref:`theory <Theory>`, the key to any CAP calculation is to 
find the optimal value of the CAP strength parameter :math:`\eta_{opt}`. The :class:`~pyopencap.analysis.EigenvalueTrajectory`
class has a function :class:`~pyopencap.analysis.EigenvalueTrajectory.find_eta_opt` for
exactly that purpose. For uncorrected trajectories, :math:`\eta_{opt}` is calculated as 

:math:`min |\eta\frac{dE}{d\eta}|`

For corrected trajectories, :math:`\eta_{opt}` is calculated as 

:math:`min |\eta\frac{dU}{d\eta}|`

The derivative is calculated numerically by means of finite differences using the gradient 
function of numpy_.

.. _numpy: https://numpy.org/doc/stable/reference/generated/numpy.gradient.html

.. code-block:: python

	uc_energy,uc_eta_opt = traj.find_eta_opt()
	corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

The presence of nonphysical stationary points can sometimes result in this function returning 
smaller values of :math:`\eta_{opt}` than desired. One can specify the `start_idx` keyword 
argument to begin the search starting from a specified index along the trajectory.

.. code-block:: python

	uc_energy,uc_eta_opt = traj.find_eta_opt(start_idx=20)
	

Visualization
--------------
It can often be helpful to visualize trajectories graphically. In addition to access to the 
class attributes `uncorrected_energies` and `corrected_energies`, we also provide some 
helper functions which process the data in useful ways for visualization. 

For instance, :func:`~pyopencap.analysis.EigenvalueTrajectory.energies_ev` returns 
the excitation energies in eV with respect to specified reference energy.

.. code-block:: python

	import matplotlib.pyplot as plt
	UC_ev = traj.energies_ev(ref_energy)
	Corr_ev = traj.energies_ev(ref_energy,corrected=True)
	plt.plot(np.real(UC_ev),np.imag(UC_ev),'ro', label='Uncorrected')
	plt.plot(np.real(Corr_ev),np.imag(Corr_ev),'ro', label='Corrected')
	plt.show()

:func:`~pyopencap.analysis.EigenvalueTrajectory.get_logarithmic_velocities` returns the value of 
:math:`\eta\frac{\partial E(\eta) }{\partial \eta}` (or 
:math:`|\eta\frac{dU}{d\eta}|\rightarrow min` if the corrected keyword argument is set to True)
for each point along the trajectory.

.. code-block:: python
	
	derivs = traj.get_logarithmic_velocities()
	plt.plot(traj.etas,derivs)
	plt.show()


References
-----------
.. [Cederbaum2002] Santra, R.; Cederbaum, L. S. Non-Hermitian Electronic Theory and Applications to Clusters. *Phys. Rep.* **2002**, 368 (1), 1–117.


.. autoclass:: pyopencap.analysis.EigenvalueTrajectory
   :members:

   .. automethod:: __init__