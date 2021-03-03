import pyopencap
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
from CAPTrajectory import CAPHamiltonian

#Change these lines to suit your system
##########################################
guess=3
eta_list = np.linspace(0,500,101)
eta_list = eta_list * 1E-5
ref_energy= -109.35498051
au2eV= 27.2113961
RASSI_FILE = "../../opencap/openmolcas/symm.rassi.h5"
OUTPUT_FILE = "../../opencap/openmolcas/symm.out"
##########################################

sys_dict = {"molecule": "molcas_rassi",
"basis_file": RASSI_FILE}

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

es_dict = {"method" : "ms-caspt2",
    "molcas_output":OUTPUT_FILE ,
        "rassi_h5":RASSI_FILE,
}

# When symmetry is turned on, data must be read in using read_data
s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,10,"openmolcas")
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()
h0 = pc.get_H()

CAPH = CAPHamiltonian(H0=h0,W=mat)
CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.all_energies),np.imag(CAPH.all_energies),'ro',label='Uncorrected Trajectory')
plt.show()

for i in range(0,CAPH.nstates):
    traj = CAPH.track_state(i,tracking="overlap")
    uc_energy,uc_eta_opt = traj.find_eta_opt()
    plt.plot(np.real(traj.uncorrected_energies),np.imag(traj.uncorrected_energies),'-ro',label='Uncorrected Trajectory')
    plt.plot(np.real(traj.corrected_energies),np.imag(traj.corrected_energies),'-bo',label='Corrected Trajectory')
    plt.legend()
    plt.show()

corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

print((uc_energy-ref_energy)*27.2114)
print((corr_energy-ref_energy)*27.2114)

plt.plot(np.real(traj.uncorrected_energies),np.imag(traj.uncorrected_energies),'-ro',label='Uncorrected Trajectory')
plt.plot(np.real(traj.corrected_energies),np.imag(traj.corrected_energies),'-bo',label='Corrected Trajectory')
plt.legend()
plt.show()
