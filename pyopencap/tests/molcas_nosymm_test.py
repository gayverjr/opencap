import pyopencap as pycap
import numpy as np
import h5py
import os
import sys

destDir="./opencap/tests/openmolcas"
sys_dict = {"molecule": "molcas_rassi",
"basis_file": destDir+"/nosymm.h5"}
molden_dict = {"molecule": "molden",
"basis_file": destDir+"/nosymm.molden"}
cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}
es_dict = {"method" : "ms-caspt2",
           "molcas_output":destDir+"/nosymm.out",
           "rassi_h5":destDir+"/nosymm.h5",
}
s1 = pycap.System(sys_dict)
s2 = pycap.System(molden_dict)
f = h5py.File(destDir+"/nosymm.h5", 'r')
nbasis = 46

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

def test_overlap():
    arr = np.array(f["AO_OVERLAP_MATRIX"])
    arr = np.reshape(arr,(nbasis,nbasis))
    s1.check_overlap_mat(arr,"openmolcas",destDir+"/nosymm.h5")

def test_tdms():
    arr = f["SFS_TRANSITION_DENSITIES"]
    pc = pycap.Projected_CAP(s1,cap_dict,10,"openmolcas")
    for i in range(0,10):
        for j in range(i,10):
            arr1 = 0.5*np.reshape(arr[i][j],(nbasis,nbasis))
            pc.add_tdms(arr1,arr1,i,j,"openmolcas",destDir+"/nosymm.h5")
            if i!=j:
                pc.add_tdms(arr1,arr1,j,i,"openmolcas",destDir+"/nosymm.h5")
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()

def test_tdm():
    arr = f["SFS_TRANSITION_DENSITIES"]
    pc = pycap.Projected_CAP(s1,cap_dict,10,"openmolcas")
    for i in range(0,10):
        for j in range(i,10):
            arr1 = np.reshape(arr[i][j],(nbasis,nbasis))
            pc.add_tdm(arr1,i,j,"openmolcas",destDir+"/nosymm.h5")
            if i!=j:
                pc.add_tdm(arr1,j,i,"openmolcas",destDir+"/nosymm.h5")
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()





