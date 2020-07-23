Tutorial
=======================
This is a tutorial to get you started using PyOpenCAP. It is an adaptation of a notebook located 
in /examples/molcas/OpenCAP.ipnyb which computes the energy and width of the pi* resonance 
of N2- using CAP-MS-CASPT2.

This specific example uses OpenMolcas as the electronic structure package. PySCF can also be used, and
the interface is very similar. More details on PySCF can be found on our PySCF page, and 
a working example can be found in /example/pyscf/h2.py.

**Preliminary: Importing the module**

In addition to PyOpenCAP, we'll also import numpy to help us process the data.

    >>> import pyopencap
    >>> import pycap
    
**Constructing the system object** 

The System object of PyOpenCAP contains the geometry and basis set information, as well
as the overlap matrix. The constructor exposed to Python takes in a Python dictionary as an argument, 
with key/value pairs which mimic the input file format of the command line version. There
are three equivalent ways of specifying the geometry and basis set: rassi_h5, molden, and 
explicit. Here, we'll use the rassi_h5 file.

    >>> sys_dict = {"molecule": "molcas_rassi","basis_file": "anion_ms.rassi.h5"}
    >>> s = pycap.System(sys_dict)
    >>> smat = s.get_overlap_mat()
    >>> np.shape(smat)
    Number of basis functions:119
    (119, 119)

**Constructing the Projected_CAP object**

The CAP matrix is computed by the "Projected_CAP" object. The constructor for the Projected_CAP 
object requires a System object, a dictionary containing the CAP parameters, 
the number of states (10 in this case), and finally the string "openmolcas", which
denotes the ordering of the atomic orbital basis set. 

    >>> cap_dict = {"cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"}
    >>> pc = pycap.Projected_CAP(s,cap_dict,10,"openmolcas")
    
**Parsing electronic structure data from file**

The "read_data" function of the Projected_CAP object can read in the effective Hamiltonian
and densities in one-shot when passed a Python dictionary with the right keywords. For now,
we'll retrieve the effective Hamiltonian and store it as h0 for later use.

    >>> es_dict = {"method" : "ms-caspt2",
           "molcas_output":"anion_ms.out",
           "rassi_h5":"anion_ms.rassi.h5"}
    >>> pc.read_data(es_dict)
    >>> h0 = pc.get_H()
    Successfully read in zeroth order Hamiltonian from file:anion_ms.out
    Successfully read in densities from file:anion_ms.rassi.h5


**Passing densities in RAM**

Alternatively, one can load in the densities one at a time using the "add_tdms or add_tdm"
functions. We load in the matrices from rassi.h5 using the h5py package, and then
pass them as numpy arrays to the Projected_CAP object. This can be particularly useful 
if we want to exclude some of the states from the Projected CAP calculation. In this 
example, the CAP matrix is made to be symmetric.

    >>> import h5py
    >>> f = h5py.File('anion_ms.rassi.h5', 'r')
    >>> dms = f["SFS_TRANSITION_DENSITIES"]
    >>> pc = pycap.Projected_CAP(s,cap_dict,10,"openmolcas")
    >>> for i in range(0,10):
    >>>     for j in range(i,10):
    >>>         dm_alpha = 0.5*np.reshape(dms[i][j],(119,119))
    >>>         dm_beta =  0.5*np.reshape(dms[i][j],(119,119))
    >>>         pc.add_tdms(dm_alpha,dm_beta,i,j,"openmolcas","anion_ms.rassi.h5")
    >>>      if i!=j:
    >>>          pc.add_tdms(dm_alpha,dm_beta,j,i,"openmolcas","anion_ms.rassi.h5")
    

Once all of the densities are loaded, the CAP matrix is computed 
using the "compute_projected_cap" function. The matrix can be retrieved using the
"get_projected_cap" function.

    >>> pc.compute_projected_cap()
    >>> W_mat=pc.get_projected_cap()

Now that we have our zeroth order Hamiltonian (stored in h0) and our CAP matrix(W_mat) in the same
wavefunction basis, we can apply our CAP as a first order perturbation by diagonalizing the full 
CAP Hamiltonian H=H0-inW. To get the best estimate of resonance position and width, we must
perform this over a range of eta values to identify the eta opt. Please see the full notebook
in /examples/molcas/OpenCAP.ipnyb or the script "cap_trajectory.py" to help get you started 
on performing the trajectory analysis. Automated tools for analyzing eigenvalue trajectories
are a subject of future work.



