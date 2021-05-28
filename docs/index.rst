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

.. _license: https://github.com/gayverjr/opencap/blob/master/LICENSE

Supported Packages
===================

* OpenMolcas_
* PySCF_
* Q-Chem_
* Psi4_

Supported Methods
=================

* EOM-CC
* FCI
* (X)-MS-CASPT2
* TDDFT
* ADC (through ADCC_)

Supported Potentials
====================

* Box
* Smooth Voronoi

Please see the :ref:`keywords <keywords>` section for more details.

Other features
==============

* Python based eigenvalue trajectory :ref:`analysis <analysis>` tools

In development
==================
* Feshbach projection approaches

.. _OpenMolcas: https://molcas.gitlab.io/OpenMolcas/sphinx/
.. _PySCF: http://pyscf.org/
.. _Psi4: http://www.psicode.org/
.. _OpenCAP: https://github.com/gayverjr/opencap
.. _pybind11: https://github.com/pybind/pybind11
.. _Q-Chem: https://www.q-chem.com/
.. _ADCC: https://adc-connect.org/

Contents
========

.. toctree::
   :maxdepth: 1

   install.rst
   getting_started.rst
   theory.rst
   interfaces.rst
   analysis/analysis.rst
   classes.rst
   keywords.rst


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`