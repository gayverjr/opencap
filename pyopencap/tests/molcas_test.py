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

def test_molcas():
    sys = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(sys,cap_dict,10)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    from pyopencap.analysis import CAPHamiltonian
    caph = CAPHamiltonian(pc=pc)
    eta_list = np.linspace(0,1500,101)
    eta_list = np.around(eta_list * 1E-5,decimals=5)
    caph.run_trajectory(eta_list)
    traj = caph.track_state(1,tracking="overlap")
    uc_energy, eta_opt = traj.find_eta_opt(start_idx=10,ref_energy=-109.35465184,units="eV")
    assert np.isclose(np.real(uc_energy),2.54772)









