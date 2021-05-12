import pyopencap
import numpy as np
from pandas import DataFrame
import h5py


FCHK_FILE = "n2.fchk"
OUTPUT_FILE = "n2.out"

sys_dict = {"molecule": "qchem_fchk",
"basis_file": FCHK_FILE}

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

es_dict = {"method" : "eom",
           "qchem_output":OUTPUT_FILE,
"qchem_fchk":FCHK_FILE}

s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,5,"qchem")
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_perturb_cap()
W=pc.get_perturb_cap()
H0 = pc.get_H()

# now its time for trajectories
from pyopencap.analysis.CAPTrajectory import CAPHamiltonian
import matplotlib.pyplot as plt
ref_energy =  -109.36195558
eta_list = np.linspace(0,5000,101)
eta_list = eta_list * 1E-5
CAPH = CAPHamiltonian(H0=H0,W=W,blah="blah")
CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.energies_ev(ref_energy=ref_energy)),np.imag(CAPH.energies_ev(ref_energy=ref_energy)),'ro',label='Uncorrected Trajectory')
plt.show()

traj = CAPH.track_state(1,tracking="overlap")
uc_energies = traj.energies_ev(ref_energy=ref_energy)
corr_energies = traj.energies_ev(ref_energy=ref_energy,corrected=True)

uc_energy, eta_opt = traj.find_eta_opt(start_idx=10)
#uc_energy = (uc_energy-ref_energy)*27.2114
corr_energy, corr_eta_opt = traj.find_eta_opt(corrected=True,start_idx=10)
#corr_energy = (corr_energy-ref_energy)*27.2114

print("Uncorrected:")
print(uc_energy)
print(eta_opt)
print("Corrected:")
print(corr_energy)
print(corr_eta_opt)

plt.plot(np.real(uc_energies),np.imag(uc_energies),'-ro',label="Uncorrected")
plt.plot(np.real(corr_energies),np.imag(corr_energies),'-bo',label="Corrected")
plt.plot(np.real(uc_energy),np.imag(uc_energy),'g*',markersize=20)
plt.plot(np.real(corr_energy),np.imag(corr_energy),'g*',markersize=20)
plt.legend()
plt.show()



