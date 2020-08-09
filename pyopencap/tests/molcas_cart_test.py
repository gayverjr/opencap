import pyopencap as pycap
import numpy as np
import h5py
import os
import sys

destDir="../opencap/tests/openmolcas"
sys_dict = {"molecule": "molcas_rassi",
"basis_file": destDir+"/xms_cart.rassi.h5"}
molden_dict = {"molecule": "molden",
"basis_file": destDir+"/xms_cart.molden"}
cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}
es_dict = {"method" : "ms-caspt2",
           "molcas_output":destDir+"/xms_cart.out",
           "rassi_h5":destDir+"/xms_cart.rassi.h5",
}
s1 = pycap.System(sys_dict)
s2 = pycap.System(molden_dict)
f = h5py.File(destDir+"/xms_cart.rassi.h5", 'r')
nbasis = 53

def test_overlap():
    arr = np.array(f["AO_OVERLAP_MATRIX"])
    arr = np.reshape(arr,(nbasis,nbasis))
    s1.check_overlap_mat(arr,"openmolcas",destDir+"/xms_cart.rassi.h5")

def test_rassi():
    pc = pycap.Projected_CAP(s1,cap_dict,3,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_molden():
    arr = np.array(f["AO_OVERLAP_MATRIX"])
    arr = np.reshape(arr,(nbasis,nbasis))
    pc = pycap.Projected_CAP(s2,cap_dict,3,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize_cap(arr,"openmolcas",destDir+"/xms_cart.rassi.h5")
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()






