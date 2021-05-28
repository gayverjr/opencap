import pyopencap
import numpy as np
import matplotlib.pyplot as plt
from pyopencap.analysis import CAPHamiltonian

ref_energy = -109.36009153
eta_list = np.linspace(0,2000,101)
eta_list = eta_list * 1E-5
RASSI_FILE = "xms.out"
OUTPUT_FILE = "xms.rassi.h5"

sys_dict = {"molecule": "molcas_rassi",
"basis_file": RASSI_FILE}

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
}

es_dict = { "package": "openmolcas",
"method" : "xms-caspt2",
    "molcas_output":OUTPUT_FILE ,
        "rassi_h5":RASSI_FILE,
}


s = pyopencap.System(sys_dict)
pc = pyopencap.CAP(s,cap_dict,10)
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_projected_cap()
W = pc.get_projected_cap()
h0 = pc.get_H()


CAPH = CAPHamiltonian(H0=h0,W=W)
CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.energies_ev(ref_energy=ref_energy)),np.imag(CAPH.energies_ev(ref_energy=ref_energy)),'ro',label='Uncorrected Trajectory')
plt.show()

for i in range(0,CAPH.nstates):
    traj = CAPH.track_state(i,tracking="overlap")
    uc_energies = traj.energies_ev(ref_energy=ref_energy)
    corr_energies = traj.energies_ev(ref_energy=ref_energy,corrected=True)
    plt.plot(np.real(uc_energies),np.imag(uc_energies),'-ro',label="Uncorrected")
    plt.plot(np.real(corr_energies),np.imag(corr_energies),'-bo',label="Corrected")
    plt.title("State: " + str(i))
    plt.legend()
    plt.show()

