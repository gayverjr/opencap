PyOpenCAP Documentation
=======================

PyOpenCAP is the Python API for OpenCAP, an open-source software aimed at extending the functionality
of modern quantum chemistry packages to describe resonances. PyOpenCAP uses the pybind11 library
to expose OpenCAP classes and methods, allowing calculations to be driven within a Python interpreter. 

PyOpenCAP is currently capable of generating the matrices required to 
perform projected complex absorbing potential calculations on metastable electronic states. 
These calculations are able to extract resonance position and width at the 
cost of a single bound-state electronic structure calculation. We currently support interfaces with 
the OpenMolcas and Pyscf quantum chemistry packages. 

To get started, please see our tutorial. For more detailed information on theory, supported
interfaces, and the python API, please see those respective sections.

*Current features*

*Upcoming features*

* automated python-based trajectory analysis tools
* interface to Psi4 quantum chemistry package


Contents:

.. toctree::
   :maxdepth: 1

   tutorial.rst
   theory.rst
   interfaces.rst
   classes.rst
   keywords.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`