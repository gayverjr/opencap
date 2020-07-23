Keywords
=======================
PyOpenCAP heavily relies on the use of Python dictionaries which contain key/value pairs
which mimic the input file format of the command line version. Here, we outline the valid 
key/value combinations. Importantly, **all key value pairs should be specified as strings**.

**System keywords**

The System object contains the basis set and geometry information, which can be obtained
in a few different ways.

+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Keyword          | Required | Default/valid values | Description                                                                                                                                                         |
+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| molecule         | yes      | molden,qchem_fchk    | Specifies which format to read the molecular geometry. If "inline" is chosen,                                                                                       |
|                  |          | rassi_h5,inline      | the "geometry" keyword is also required.                                                                                                                            |
+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| geometry         | no       | See below            | Specifies the geometry in an inline format described below. Required when the                                                                                       |
|                  |          |                      | "molecule" field is set to "inline".                                                                                                                                |
+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| basis_file       | yes      | path to basis file   | Specifies the path to the basis file. When "molecule" is set to "molden","rassi_h5", or "qchem_fchk",                                                               |
|                  |          |                      | this field should be set to a path to a file of the specified type. When "molecule" is set to                                                                       |
|                  |          |                      | "inline", this field should be set to a path to a basis set file formatted in "Psi4" style.                                                                         |
+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| cart_bf          | no       | ""                   | Controls the use of pure or Cartesian angular forms of GTOs. The letters corresponding to the angular momenta listed in this field will be expanded in cartesians,  |
|                  |          |                      | those not listed will be expanded in pure GTOs. For example, "df" means d and f-type                                                                                |
|                  |          |                      | functions will be cartesian, and all others will be pure harmonic. This keyword is only active                                                                      |
|                  |          |                      | when "molecule" is set to "inline".                                                                                                                                 |
+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| bohr_coordinates | no       | "false"              | Set to true when coordinates specified in "geometry" keyword are in bohr units. This keyword is only active when "molecule" is set to "inline".                     |
+------------------+----------+----------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------+

When specifying the geometry inline, use the following format:

atom1 x-coordinate y-coordinate z-coordinate 

atom2 x-coordinate y-coordinate z-coordinate ...

Ghost centers with zero nuclear charge can be specified using the symbol "X".

Example:

.. code-block:: python

    sys_dict = {"geometry":    '''N  0  0   1.039
                              N  0  0   -1.039
                              X   0  0   0.0''',
            		"molecule" : "read",
            		"basis_file":"path/to/basis.bas",
            		"cart_bf":"d",
            		"bohr_coordinates:": "true"}

**CAP keywords**

PyOpenCAP supports Voronoi and Box-type absorbing potentials. We also allow some customization
of the numerical grid used for integration. Please see https://github.com/dftlibs/numgrid for
more details on the radial_precision and angular_points keywords.

+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| Keyword          | Required | Default/valid values | Description                                                                                                        |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| cap_type         | yes      | box or voronoi       | Type of absorbing potential.                                                                                       |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| cap_x            | no       | N/A                  | Onset of CAP in x-direction. Specify in bohr units. Only active when "cap_type' is set to "box".                   |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| cap_y            | no       | N/A                  | Onset of CAP in y-direction. Specify in bohr units. Only active when "cap_type' is set to "box".                   |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| cap_y            | no       | N/A                  | Onset of CAP in z-direction. Specify in bohr units. Only active when "cap_type' is set to "box".                   |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| r_cut            | no       | N/A                  | Cutoff radius for Voronoi CAP. Specify in bohr units. Only active when "cap_type' is set to "voronoi".             |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| radial_precision | no       | 14                   | Radial precision for numerical integration grid. A precision of 1x10^(-N), where N is the value specified is used. |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+
| angular_points   | no       | 590                  | Number of angular points used for the grid. See https://github.com/dftlibs/numgrid for allowed numbers of points.  |
+------------------+----------+----------------------+--------------------------------------------------------------------------------------------------------------------+

Example:

.. code-block:: python

    cap_dict = {"cap_type": "box",
            	"cap_x":"2.76",
            	"cap_y":"2.76",
            	"cap_z":"4.88",
            	"Radial_precision": "14",
            	"angular_points": "110"}

**Electronic structure keywords**

The "read_data" function of the Projected_CAP class is able to parse the zeroth order Hamiltonian
and load the densities when supplied with an appropriate formatted dictionary. All keywords
must be specified to use this function. Currently, this is only supported for calculations
using the OpenMolcas interface.

+---------------+------------------------------------------------------------------------------------------------------+
| Keyword       | Description                                                                                          |
+---------------+------------------------------------------------------------------------------------------------------+
| method        | Electronic structure method used in the calculation. Valid options are "MS-CASPT2" and "XMS-CASPT2". |
+---------------+------------------------------------------------------------------------------------------------------+
| molcas_output | Path to OpenMolcas output file.                                                                      |
+---------------+------------------------------------------------------------------------------------------------------+
| rassi_h5      | Path to OpenMolcas rassi.h5 file.                                                                    |
+---------------+------------------------------------------------------------------------------------------------------+

Example:

.. code-block:: python
	
    es_dict = {"method" : "ms-caspt2",
           "molcas_output":"path/to/output.out",
           "rassi_h5":"path/to/rassi.h5"}
    pc.read_data(es_dict)



