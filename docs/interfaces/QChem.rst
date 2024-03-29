QChem
=======
PyOpenCAP supports an interface with the Q-Chem_ quantum chemistry package. 

.. _Q-Chem: https://www.q-chem.com/

Importing data
---------------

**System object**

The geometry and basis set can be imported into a :class:`~pyopencap.System` 
object using .fchk files.

.. code-block:: python

	import pyopencap
	sys_dict = {"molecule": "qchem_fchk","basis_file": "path/to/qc.fchk"}
	my_system = pyopencap.System(sys_dict)

**CAP object**
Densities can be read in from .fchk files, and the zeroth order Hamiltonian can be read from 
Q-Chem output files for EOM-CC calculations. To export the full densities to .fchk, GUI=2 
must be set in the $rem card, and PROJ_CAP=3 must be set in the 
$complex_ccman card. See the Q-Chem manual_ for more details.

.. _manual: https://manual.q-chem.com/latest/sub_cc_cs_cap.html

The following snippet can be used to read the data from a Q-Chem output and properly formatted 
.fchk file, and calculate the CAP matrix:

.. code-block:: python

    cap_dict = {"cap_type": "box",
            	"cap_x":"2.76",
            	"cap_y":"2.76",
            	"cap_z":"4.88",
            	"Radial_precision": "14",
            	"angular_points": "110"}
    nstates = 10
    pc = pyopencap.CAP(my_system,cap_dict,nstates)
    # read in densities
	es_dict = {"method" : "eom",
	           "package": "qchem",
           "qchem_output":"path/to/output.out",
           "qchem_fchk":"path/to/qc.fchk"}
    pc.read_data(es_dict)
    # save the zeroth order Hamiltonian for later use
    h0 = pc.get_H()
    pc.compute_projected_cap()
    W_mat=pc.get_projected_cap()
    
Generate and analyze eigenvalue trajectories
-----------------------------------------------
H0 and W, or the :class:`~pyopencap.CAP` object can be used to construct a :class:`~pyopencap.analysis.CAPHamiltonian` object. 

.. code-block:: python

	from pyopencap.analysis import CAPHamiltonian
	CAPH = CAPHamiltonian(H0=H0,W=W_mat)
	# equivalently
	CAPH = CAPHamiltonian(pc=pc)

Additionally, Q-Chem (starting from version 5.4) natively implements Projected CAP-EOM-CC and 
Projected CAP-ADC methods, and prints the necessary matrices to the output. PyOpenCAP can parse 
these output files to generate :class:`~pyopencap.analysis.CAPHamiltonian` objects.

.. code-block:: python

	from pyopencap.analysis import CAPHamiltonian
	CAPH = CAPHamiltonian(output="proj-eomcc.out",irrep="B2g")
	CAPH = CAPHamiltonian(output="proj-adc.out",onset="3000")

See the :ref:`analysis <analysis>` section for more details.
	
