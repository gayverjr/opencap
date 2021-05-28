# Analyze Q-Chem output and checkpoint file
import pyopencap
from pyopencap.analysis import CAPHamiltonian
import matplotlib.pyplot as plt
import numpy as np

ref_energy = -109.36195558

sys_dict = {"molecule":  "qchem_fchk",
            "basis_file": "qc_inp.fchk"
}

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88"
}

es_dict = { "package": "qchem",
            "method" : "eom",
           "qchem_output":"qc_inp.out",
           "qchem_fchk":"qc_inp.fchk",
}

s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,5)
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_projected_cap()
W = pc.get_projected_cap()
H0 = pc.get_H()

CAPH = CAPHamiltonian(H0=H0,W=W)
eta_list = np.linspace(0,5000,101) * 1E-5
CAPH.run_trajectory(eta_list)


for i in range(0,CAPH.nstates):
    traj = CAPH.track_state(i,tracking="overlap")
    uc_energies = traj.energies_ev(ref_energy=ref_energy)
    corr_energies = traj.energies_ev(ref_energy=ref_energy,corrected=True)
    plt.plot(np.real(uc_energies),np.imag(uc_energies),'-ro',label="Uncorrected")
    plt.plot(np.real(corr_energies),np.imag(corr_energies),'-bo',label="Corrected")
    plt.title("State: " + str(i))
    plt.legend()
    plt.show()


traj = CAPH.track_state(1,tracking="overlap")
uc_energies = traj.energies_ev(ref_energy=ref_energy)
corr_energies = traj.energies_ev(ref_energy=ref_energy,corrected=True)
uc_energy_au,uc_eta = traj.find_eta_opt(ref_energy=ref_energy,units="au")
corr_energy_au,corr_eta = traj.find_eta_opt(ref_energy=ref_energy,corrected=True,units="au")
uc_energy_ev,uc_eta = traj.find_eta_opt(ref_energy=ref_energy,units="eV")
corr_energy_ev,corr_eta = traj.find_eta_opt(ref_energy=ref_energy,corrected=True,units="eV")
plt.plot(np.real(uc_energies),np.imag(uc_energies),'-ro',label="Uncorrected")
plt.plot(np.real(corr_energies),np.imag(corr_energies),'-bo',label="Corrected")
plt.plot(np.real(uc_energy),np.imag(uc_energy),'g*',markersize=15)
plt.plot(np.real(corr_energy),np.imag(corr_energy),'g*',markersize=15)
plt.legend()
plt.show()



