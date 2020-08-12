pyopencap.CAP
========================
The :class:`~pyopencap.CAP` class is used to compute the CAP matrix first in AO basis, and then in 
wave function basis using the one-particle densities which are passed in. It is also capable of 
parsing OpenMolcas output files to obtain the zeroth order Hamiltonian and return it to the user.

.. currentmodule:: pyopencap

.. autoclass:: CAP

   .. automethod:: __init__
   .. automethod:: add_tdm
   .. automethod:: add_tdms
   .. automethod:: compute_ao_cap
   .. automethod:: compute_perturb_cap
   .. automethod:: get_H
   .. automethod:: get_ao_cap
   .. automethod:: get_perturb_cap  
   .. automethod:: read_data 
   .. automethod:: renormalize
   .. automethod:: renormalize_cap