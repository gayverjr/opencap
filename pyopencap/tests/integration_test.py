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
    SOFTWARE.
    '''

import pyopencap
from pyopencap.analysis import CAPHamiltonian
import numpy as np
import os
import sys
import pytest

try:
    import numgrid
except:
    pass

def box_cap(x,y,z,w):
    cap_values = []
    cap_x = 3.00
    cap_y = 3.00
    cap_z = 3.00
    for i in range(0,len(x)):
        result = 0
        if np.abs(x[i])>cap_x:
            result += (np.abs(x[i])-cap_x) * (np.abs(x[i])-cap_x)
        if np.abs(y[i])>cap_y:
            result += (np.abs(y[i])-cap_y) * (np.abs(y[i])-cap_y)
        if np.abs(z[i])>cap_z:
            result += (np.abs(z[i])-cap_z) * (np.abs(z[i])-cap_z)
        result = w[i]*result
        cap_values.append(result)
    return cap_values

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

cap_dict = {
            "cap_type": "box",
            "cap_x":"3.00",
            "cap_y":"3.00",
            "cap_z":"3.00",
}

cap_dict_numerical = {
            "cap_type": "box",
            "cap_x":"3.00",
            "cap_y":"3.00",
            "cap_z":"3.00",
            "do_numerical": "true"
}

ref_energy = -109.36195558

radial_precision = 1.0e-12
min_num_angular_points = 590
max_num_angular_points = 590
proton_charges = [7, 7, 1]
center_coordinates_bohr = [(0.0, 0.0, 1.03699997), (0.0, 0.0, -1.03699997), (0.0, 0.0, 0.0)]

def get_corrected_energy(pc):
    CAPH = CAPHamiltonian(pc=pc)
    eta_list = np.linspace(0,5000,101)
    eta_list = np.around(eta_list * 1E-5,decimals=5)
    CAPH.run_trajectory(eta_list)
    traj = CAPH.track_state(1,tracking="energy")
    corr_energy, corr_eta_opt = traj.find_eta_opt(corrected=True,start_idx=10,end_idx=-1,ref_energy=ref_energy,units="eV")
    return corr_energy

def test_custom_cap():
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,{"cap_type": "custom"},5,box_cap)
    pc.read_data(es_dict)
    pc.compute_projected_cap()
    e1  = get_corrected_energy(pc)
    pc.compute_ao_cap(cap_dict)
    pc.compute_projected_cap()
    e2  = get_corrected_energy(pc)
    assert np.isclose(e1,e2)

def test_numerical_integration():
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,cap_dict,5)
    pc.read_data(es_dict)
    pc.compute_projected_cap()
    e1  = get_corrected_energy(pc)
    pc.compute_ao_cap(cap_dict_numerical)
    e2  = get_corrected_energy(pc)
    assert np.isclose(e1,e2)

@pytest.mark.skipif('numgrid' not in sys.modules,reason="requires the Numgrid library")
def test_custom_grid():
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,cap_dict,5)
    pc.read_data(es_dict)
    pc.compute_projected_cap()
    e1  = get_corrected_energy(pc)
    pc = pyopencap.CAP(s,cap_dict,5)
    pc.read_data(es_dict)    
    import numgrid
    for center_index in range(len(center_coordinates_bohr)):
        coordinates, w = numgrid.atom_grid_bse(
            "aug-cc-pvtz",
            1.0e-12,
            590,
            590,
            proton_charges,
            center_index,
            center_coordinates_bohr,
            hardness=3,
        )
        coordinates = np.array(coordinates)
        x = np.array(coordinates[:,0])
        y = np.array(coordinates[:,1])
        z = np.array(coordinates[:,2])
        w = np.array(w)
        pc.compute_cap_on_grid(x,y,z,w)
    pc.compute_projected_cap()
    e2  = get_corrected_energy(pc)
    assert np.isclose(e1,e2)