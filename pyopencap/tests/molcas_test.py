'''Copyright (c) 2021 James Gayvert
    
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


# Benchmark for OpenMolcas

import pyopencap
import os
from pyopencap.analysis import CAPHamiltonian
import numpy as np

dest_dir = "../examples/analysis/N2/ref_outputs"
RASSI_FILE = os.path.join(dest_dir,"xms.rassi.h5")
OUTPUT_FILE = os.path.join(dest_dir,"xms.out")
h0_file = "../examples/opencap/heff.in"

sys_dict = {
    "molecule": "molcas_rassi",
    "basis_file":RASSI_FILE
}

cap_dict ={
            "cap_type": "voronoi",
            "r_cut":"3.00",
            "Radial_precision": "14",
            "angular_points": "110"
}

es_dict ={
    "package":"openmolcas",
    "method":"xms-caspt2",
    "molcas_output":OUTPUT_FILE ,
    "rassi_h5":RASSI_FILE
}

es_dict2 ={
    "package":"openmolcas",
    "method":"ms-caspt2",
    "h0_file":h0_file,
    "rassi_h5":RASSI_FILE
}

def test_molcas():
    sys = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(sys,cap_dict,10)
    pc.read_data(es_dict)
    pc.compute_ao_cap(cap_dict)
    pc.compute_projected_cap()
    CAPH = CAPHamiltonian(pc=pc)
    

def test_heff():
    sys = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(sys,cap_dict,10)
    pc.read_data(es_dict2)
    h0 = pc.get_H()
    print(h0[0][0])
    assert np.isclose(h0[0][0],-109.312105)









