import pyopencap
import numpy as np
import h5py
import os
import sys

destDir="../opencap/tests/openmolcas"
sys_dict = {"molecule": "molcas_rassi",
"basis_file": destDir+"/symm.rassi.h5"}
molden_dict = {"molecule": "molden",
"basis_file": destDir+"/symm.molden"}
cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}
es_dict = {"method" : "ms-caspt2",
           "molcas_output":destDir+"/symm.out",
           "rassi_h5":destDir+"/symm.rassi.h5",
}
s1 = pyopencap.System(sys_dict)
s2 = pyopencap.System(molden_dict)
f = h5py.File(destDir+"/symm.rassi.h5", 'r')

def test_rassi():
    pc = pyopencap.CAP(s1,cap_dict,3,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()
    h0 = pc.get_H()

def test_molden():
    pc = pyopencap.CAP(s2,cap_dict,3,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()
    h0 = pc.get_H()







