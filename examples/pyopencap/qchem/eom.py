# Analyze Q-Chem output and checkpoint file
import pyopencap
from pyopencap.analysis import CAPHamiltonian
import matplotlib.pyplot as plt

sys_dict = {"molecule":  "qchem_fchk",
            "basis_file": "n2.fchk"
}

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88"
}

es_dict = { "package": "qchem",
            "method" : "eom",
           "qchem_output":"n2.out",
           "qchem_fchk":"n2.fchk",
}

s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,5)
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_projected_cap()
W = pc.get_projected_cap()
H0 = pc.get_H()

CAPH = CAPHamiltonian(H0=H0,W=W)
eta_list = np.linspace(0,2000,101) * 1E-5
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

