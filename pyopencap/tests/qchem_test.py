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

import pyopencap
import os 
import numpy as np


cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

dest_dir = "../examples/analysis/N2/ref_outputs"
FCHK_FILE = os.path.join(dest_dir,"qc_inp.fchk")
OUTPUT_FILE = os.path.join(dest_dir,"qc_inp.out")
h0_file = "../examples/opencap/heff_diag.in"

sys_dict = {
    "molecule": "qchem_fchk",
    "basis_file":FCHK_FILE
}

es_dict ={
    "package":"qchem",
    "method":"eom",
    "qchem_output":OUTPUT_FILE,
    "qchem_fchk":FCHK_FILE,
}

es_dict2 ={
    "package":"qchem",
    "method":"eom",
    "h0_file":h0_file,
    "qchem_fchk":FCHK_FILE,
}

def test_qchem():
    sys = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(sys,cap_dict,5)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    h0 = pc.get_H()
    from pyopencap.analysis import CAPHamiltonian
    caph = CAPHamiltonian(H0=h0,W=W)
    eta_list = np.linspace(0,5000,101)
    eta_list = np.around(eta_list * 1E-5,decimals=5)
    caph.run_trajectory(eta_list)
    traj = caph.track_state(1,tracking="overlap")
    corr_energy, eta_opt = traj.find_eta_opt(start_idx=10,ref_energy=-109.36195558,corrected=True,units="eV")
    assert np.isclose(np.real(corr_energy),2.59667)

def test_heff():
    sys = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(sys,cap_dict,5)
    pc.read_data(es_dict2)
    h0 = pc.get_H()
    assert np.isclose(h0[0][0],-109.313539)




