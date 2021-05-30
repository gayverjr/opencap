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
cap_dict =
{
            "cap_type": "voronoi",
            "r_cut":"3.00",
            "Radial_precision": "14",
            "angular_points": "110"
}

def test_ms_caspt2():
    sys_dict = {}
    es_dict = {"method" : "ms-caspt2",
"package": "openmolcas",
"molcas_output":destDir+"/ms_spherical.out",
"rassi_h5":destDir+"/ms_spherical.rassi.h5"}
    pc = pyopencap.CAP(s1,cap_dict,16)
    pc.read_data(es_dict2)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    h0 = pc.get_H()

def test_molden():
    pc = pyopencap.CAP(s2,cap_dict,3)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    h0 = pc.get_H()
    assert 

def test_pc_nevpt2():
    pc = pyopencap.CAP(s3,cap_dict,10)
    pc.read_data(es_dict3)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    h0 = pc.get_H()
    assert np.isclose(h0[0][0],-109.3332404495)

def test_sc_nevpt2():
    pc = pyopencap.CAP(s3,cap_dict,10)
    pc.read_data(es_dict4)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    h0 = pc.get_H()
    assert np.isclose(h0[0][0],-109.3327110632)




