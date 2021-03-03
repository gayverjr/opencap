import pyopencap
import numpy as np
from pandas import DataFrame
import h5py
import matplotlib.pyplot as plt
from CAPTrajectory import CAPHamiltonian

#Change these lines to suit your system
##########################################
ref_energy = -113.14485360
guess = 2.2
eta_list = np.linspace(0,5000,501)
eta_list = eta_list * 1E-5
FCHK_FILE = "../../opencap/qchem/CO.fchk"
OUTPUT_FILE = "../../opencap/qchem/CO.out"
##########################################

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

es_dict = {"method" : "eomea",
           "qchem_output":OUTPUT_FILE,
"qchem_fchk":FCHK_FILE}

# Method 1: Read data from output
s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,10,"qchem")
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()
h0 = pc.get_H()


CAPH = CAPHamiltonian(H0=h0,W=mat)
CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.energies_ev(ref_energy=ref_energy)),np.imag(CAPH.energies_ev(ref_energy=ref_energy)),'ro',label='Uncorrected Trajectory')
plt.show()

traj = CAPH.track_state(4,tracking="overlap")
uc_energy,uc_eta_opt = traj.find_eta_opt()
corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

print((uc_energy-ref_energy)*27.2114)
print((corr_energy-ref_energy)*27.2114)

plt.plot(np.real(traj.energies_ev(ref_energy=ref_energy)),np.imag(traj.energies_ev(ref_energy=ref_energy)),'-ro',label='Uncorrected Trajectory')
plt.plot(np.real(traj.energies_ev(ref_energy=ref_energy,corrected=True)),np.imag(traj.energies_ev(ref_energy=ref_energy,corrected=True)),'-bo',label='Corrected Trajectory')
plt.legend()
plt.show()
