'''Copyright (c) 2022 James Gayvert
    
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
import os
from pyopencap.analysis import CAPHamiltonian, colparser
import numpy as np

dest_dir = "../examples/pyopencap/columbus"
data_dir = 'data_files'
cap_dict = {"cap_type": "voronoi","r_cut": "3.00"}
sys_dict = {"molecule": "molden", "basis_file": os.path.join(dest_dir,data_dir,'molden_mo_mc.sp')}
parser = colparser(os.path.join(dest_dir,data_dir,'molden_mo_mc.sp'),os.path.join(dest_dir,data_dir,'tranls'))
sys = pyopencap.System(sys_dict)
pc = pyopencap.CAP(sys, cap_dict, 4)
smat = sys.get_overlap_mat('molden')


def test_h0():
    # default is pople
    H0 = parser.get_H0(filename=os.path.join(dest_dir,data_dir,'ciudgsm'))
    assert H0[0][0] == -109.331057316404
    H0 = parser.get_H0(correction_type='eci+dv2',filename=os.path.join(dest_dir,data_dir,'ciudgsm'))
    assert H0[0][0] == -109.334157554115

def test_absolute():
    sdm = parser.sdm_ao(1,filename=os.path.join(dest_dir,data_dir,'cid1trfl.FROMdrt1.state1TOdrt1.state1.iwfmt'))
    assert np.round(np.trace(smat@sdm),decimals=1) == 15.0

def test_relative():
    sdm = parser.sdm_ao(1,DRTn=1,data_dir=os.path.join(dest_dir,data_dir))
    assert np.round(np.trace(smat@sdm),decimals=1) == 15.0
    tdm = parser.tdm_ao(1,1,drtFrom=1,drtTo=1,data_dir=os.path.join(dest_dir,data_dir))
    assert np.round(np.trace(smat@tdm),decimals=1) == 11.0

def test_mo_summary():
    output = parser.mo_summary()
    assert 'Total number of basis functions: 119' in output
    assert 'NFCPSY: [1, 0, 0, 0, 1, 0, 0, 0] ' in output

def test_columbus():
    nstates = 4
    for i in range(0,nstates):
        for j in range(i,nstates):
            if i==j:
                dm1_ao = parser.sdm_ao(i+1,data_dir=os.path.join(dest_dir,data_dir))
                pc.add_tdm(dm1_ao,i,j,'molden')
            else:
                dm1_ao = parser.tdm_ao(i+1, j+1,data_dir=os.path.join(dest_dir,data_dir))
                pc.add_tdm(dm1_ao,i,j,'molden')
                pc.add_tdm(dm1_ao.conj().T,j,i,'molden')
    pc.compute_projected_cap()
    W=pc.get_projected_cap()
    assert np.sum(W) != 0.0
