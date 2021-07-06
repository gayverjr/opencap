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

cap_dict = {
    "cap_type": "voronoi",
        "r_cut":"3.00",
            "radial_precision":"14",
            "angular_points":"590"
}

es_dict = { "package": "qchem",
    "method" : "eom",
        "qchem_output":"../../analysis/N2/ref_outputs/qc_inp.out",
            "qchem_fchk":"../../analysis/N2/ref_outputs/qc_inp.fchk",
}

s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,5)
pc.read_data(es_dict)
pc.compute_projected_cap()
W = pc.get_projected_cap()
H0 = pc.get_H()

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



