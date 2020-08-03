import pyopencap as pycap
import numpy as np
import h5py
import os
import sys

destDir=sys.path[0]+"/openmolcas"
sys_dict = {"molecule": "molcas_rassi",
"basis_file": destDir+"/symm.h5"}
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
           "rassi_h5":destDir+"/symm.h5",
}
s1 = pycap.System(sys_dict)
s2 = pycap.System(molden_dict)
f = h5py.File(destDir+"/symm.h5", 'r')

def test_rassi():
    pc = pycap.Projected_CAP(s1,cap_dict,10,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_molden():
    pc = pycap.Projected_CAP(s2,cap_dict,10,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()







