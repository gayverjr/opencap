PyOpenCAP Documentation
=======================
PyOpenCAP is the Python API for OpenCAP_, an open-source software aimed at extending the functionality
of quantum chemistry packages to describe resonances. PyOpenCAP uses the pybind11_ library
to expose C++ classes and methods, allowing calculations to be driven within a Python interpreter. 

PyOpenCAP is currently capable of processing quantum chemistry data in order to
perform 'projected' complex absorbing potential calculations on metastable electronic states. 
These calculations are able to extract resonance position and width at the 
cost of a single bound-state electronic structure calculation.

To get started, please see our :ref:`Getting Started <getting_started>` page.

If you have questions or need support, please open an issue on GitHub, or contact us directly at gayverjr@bu.edu.

PyOpenCAP is released under the MIT license_.

.. _license: https://github.com/gayverjr/opencap/blob/main/LICENSE

Supported Packages
===================

* OpenMolcas_
* PySCF_
* Q-Chem_
* Psi4_
* Columbus_

Supported Methods
=================

* EOM-CC
* FCI
* MS-CASPT2 (and extended variants)
* TDDFT
* ADC (through ADCC_)
* MR-CI family of methods

Supported Potentials
====================

* Box (analytical integrals are now available!)
* Smooth Voronoi
* :ref:`Custom <Custom>` CAP functions

Please see the :ref:`keywords <keywords>` section for more details.

Other features
==============

* Python based eigenvalue trajectory :ref:`analysis <analysis>` tools
* :ref:`Custom <custom>` numerical integration grids

In development
==================
* Feshbach projection approaches
* Tools for dyanmics on complex potential energy surfaces

.. _OpenMolcas: https://molcas.gitlab.io/OpenMolcas/sphinx/
.. _PySCF: http://pyscf.org/
.. _Psi4: http://www.psicode.org/
.. _OpenCAP: https://github.com/gayverjr/opencap
.. _pybind11: https://github.com/pybind/pybind11
.. _Q-Chem: https://www.q-chem.com/
.. _ADCC: https://adc-connect.org/
.. _repository: https://github.com/tsommerfeld/L2-methods_for_resonances
.. _website: https://people.sc.fsu.edu/~jburkardt/
.. _MolSSI: https://molssi.org/
.. _Columbus: https://www.univie.ac.at/columbus/


Acknowledgements
=================

* We would like to give special thanks to Prof. Thomas Sommerfeld for generously providing us with a 
  C++ implementation of analytical box CAP integrals. Please checkout his repository_ on GitHub for many 
  implementations of L^2 methods for resonances!
* We would like to thank Prof. John Burkardt for his C++ implementation of the Incomplete Gamma function. Many of his 
  wonderful open source codes can be found at his website_.
* We would like to thank the Molecular Sciences Software Institute for funding this project. Please visit the MolSSI_ 
  website for their educational resources, fellowship opportunities, and software projects!
* We would like to thank Mushir Thodika from Temple University for his help in developing the interface with Columbus.

Contents
========

.. toctree::
   :maxdepth: 1

   install.rst
   getting_started.rst
   theory.rst
   interfaces.rst
   analysis/analysis.rst
   custom_cap_and_grid.rst
   classes.rst
   keywords.rst


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`