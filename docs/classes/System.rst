pyopencap.System
================
The System class is used to store the molecular geometry and the basis set. Upon construction,
it automatically computes the overlap matrix which can be accessed and used to verify the 
the ordering of the atomic orbital basis set.

.. currentmodule:: pyopencap

.. autoclass:: System

   
   .. automethod:: __init__
   .. automethod:: get_overlap_mat
   .. automethod:: check_overlap_mat

