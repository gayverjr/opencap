import pycap
import numpy as np
from pandas import DataFrame
import h5py

sys_dict = {"geometry":    '''N  0  0   1.039
                             N  0  0   -1.039
                            Gh 0  0   0.0''',
            "basis_file":"test_bas.bas",
            "bohr_coordinates": "true",
            "cart_bf": ""}
cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}
es_dict = {"method" : "ms-caspt2",
           "package": "openmolcas",
           "molcas_output":"anion_reference.out",
           "rassi_h5":"test.rassi.h5",
           "nstates": "10"
}

f = h5py.File('test.rassi.h5', 'r')
arr = f["SFS_TRANSITION_DENSITIES"]

nstates = 10
f = open("heff.in","r")
lines = f.readlines()
idx = 1
all_values = []
while idx<len(lines):
    vals = lines[idx].split()
    for val in vals:
        all_values.append(float(val))
    idx+=1
all_values = np.reshape(all_values,(nstates,nstates))

s = pycap.System(sys_dict)
pc = pycap.Projected_CAP(s,10,"openmolcas")
pc.read_data(es_dict)
'''
pc.set_h0(all_values)
for i in range(0,10):
    for j in range(i,10):
        arr1 = 0.5*np.reshape(arr[i][j],(119,119))
        pc.set_tdms(arr1,arr1,i,j)
        if i!=j:
            pc.set_tdms(arr1,arr1,j,i)
'''
pc.set_cap_params(cap_dict)
pc.run()
mat=pc.get_CAP_mat()
print(DataFrame(mat).to_string(index=False, header=False))

