import pyopencap
import numpy as np
import h5py
import os
import sys

destDir="../opencap/tests/openmolcas"
sys_dict = {"molecule": "molcas_rassi",
"basis_file": destDir+"/ms_spherical.rassi.h5"}
molden_dict = {"molecule": "molden",
"basis_file": destDir+"/ms_spherical.molden"}
cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}
es_dict = {"method" : "ms-caspt2",
           "molcas_output":destDir+"/ms_spherical.out",
           "rassi_h5":destDir+"/ms_spherical.rassi.h5",
}
s1 = pyopencap.System(sys_dict)
s2 = pyopencap.System(molden_dict)
f = h5py.File(destDir+"/ms_spherical.rassi.h5", 'r')
nbasis = 43

def test_rassi():
    pc = pyopencap.CAP(s1,cap_dict,16,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()
    h0 = pc.get_H()

def test_molden():
    pc = pyopencap.CAP(s2,cap_dict,16,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()
    h0 = pc.get_H()

def test_overlap():
    arr = np.array(f["AO_OVERLAP_MATRIX"])
    arr = np.reshape(arr,(nbasis,nbasis))
    s1.check_overlap_mat(arr,"openmolcas",destDir+"/ms_spherical.rassi.h5")

def test_tdms():
    arr = f["SFS_TRANSITION_DENSITIES"]
    pc = pyopencap.CAP(s1,cap_dict,16,"openmolcas")
    for i in range(0,16):
        for j in range(i,16):
            arr1 = 0.5*np.reshape(arr[i][j],(nbasis,nbasis))
            pc.add_tdms(arr1,arr1,i,j,"openmolcas",destDir+"/ms_spherical.rassi.h5")
            if i!=j:
                pc.add_tdms(arr1,arr1,j,i,"openmolcas",destDir+"/ms_spherical.rassi.h5")
    pc.compute_ao_cap()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()

def test_tdm():
    arr = f["SFS_TRANSITION_DENSITIES"]
    pc = pyopencap.CAP(s1,cap_dict,16,"openmolcas")
    for i in range(0,16):
        for j in range(i,16):
            arr1 = np.reshape(arr[i][j],(nbasis,nbasis))
            pc.add_tdm(arr1,i,j,"openmolcas",destDir+"/ms_spherical.rassi.h5")
            if i!=j:
                pc.add_tdm(arr1,j,i,"openmolcas",destDir+"/ms_spherical.rassi.h5")
    pc.compute_ao_cap()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()





