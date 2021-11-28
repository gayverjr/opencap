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

# Analyze Q-Chem output and checkpoint file
import pyopencap
from pyopencap.analysis import CAPHamiltonian
import numpy as np

sys_dict = {"molecule":  "qchem_fchk",
    "basis_file": "../../analysis/N2/ref_outputs/qc_inp.fchk"
}

es_dict = { "package": "qchem",
    "method" : "eom",
        "qchem_output":"../../analysis/N2/ref_outputs/qc_inp.out",
            "qchem_fchk":"../../analysis/N2/ref_outputs/qc_inp.fchk",
}

cap_dict = {
    "cap_type": "custom"
}

# this defines a box CAP of with cutoffs of 3 bohr in each coordinate
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

radial_precision = 1.0e-12
min_num_angular_points = 590
max_num_angular_points = 590
proton_charges = [7, 7, 1]
center_coordinates_bohr = [(0.0, 0.0, 1.03699997), (0.0, 0.0, -1.03699997), (0.0, 0.0, 0.0)]

s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,5,box_cap)
pc.read_data(es_dict)

# one can use the numgrid python API (https://github.com/dftlibs/numgrid) to allocate the grid
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
W = pc.get_projected_cap()
H0 = pc.get_H()

# now running analysis
CAPH = CAPHamiltonian(pc=pc)
eta_list = np.linspace(0,5000,201)
eta_list = np.around(eta_list * 1E-5,decimals=5)
CAPH.run_trajectory(eta_list)
traj = CAPH.track_state(1,tracking="energy")
ref_energy = -109.36195558
# Find optimal value of eta
uc_energy, eta_opt = traj.find_eta_opt(start_idx=10,ref_energy=ref_energy,units="eV")
# start_idx and end_idx for search use python slice notation (i.e. [start_idx:end_idx]).
corr_energy, corr_eta_opt = traj.find_eta_opt(corrected=True,start_idx=10,end_idx=-1,ref_energy=ref_energy,units="eV")
uc_energy_au = traj.get_energy(eta_opt,units="au")
corr_energy_au = traj.get_energy(eta_opt,units="au",corrected=True)
print("Uncorrected:")
print(uc_energy)
print(uc_energy_au)
print(eta_opt)
print("Corrected:")
print(corr_energy)
print(corr_energy_au)
print(corr_eta_opt)

import matplotlib.pyplot as plt
plt.plot(np.real(CAPH.energies_ev(ref_energy)),np.imag(CAPH.energies_ev(ref_energy)),'ro')
plt.show()
plt.plot(np.real(traj.energies_ev(ref_energy)),np.imag(traj.energies_ev(ref_energy)),'-ro')
plt.plot(np.real(traj.energies_ev(ref_energy,corrected=True)),np.imag(traj.energies_ev(ref_energy,corrected=True)),'-bo')
plt.plot(np.real(corr_energy),np.imag(corr_energy),"g*",markersize=20)
plt.plot(np.real(uc_energy),np.imag(uc_energy),"g*",markersize=20)
plt.show()



