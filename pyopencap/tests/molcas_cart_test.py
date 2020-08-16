'''Copyright (c) 2020 James Gayvert
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.'''

import pyopencap
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
s1 = pyopencap.System(sys_dict)
s2 = pyopencap.System(molden_dict)
f = h5py.File(destDir+"/xms_cart.rassi.h5", 'r')
nbasis = 53

def test_overlap():
    arr = np.array(f["AO_OVERLAP_MATRIX"])
    arr = np.reshape(arr,(nbasis,nbasis))
    s1.check_overlap_mat(arr,"openmolcas",destDir+"/xms_cart.rassi.h5")

def test_rassi():
    pc = pyopencap.CAP(s1,cap_dict,3,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize()
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()
    h0 = pc.get_H()

def test_molden():
    arr = np.array(f["AO_OVERLAP_MATRIX"])
    arr = np.reshape(arr,(nbasis,nbasis))
    pc = pyopencap.CAP(s2,cap_dict,3,"openmolcas")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize_cap(arr,"openmolcas",destDir+"/xms_cart.rassi.h5")
    pc.compute_perturb_cap()
    mat=pc.get_perturb_cap()
    h0 = pc.get_H()






