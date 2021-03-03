import pyopencap
import numpy as np
from pandas import DataFrame
import h5py
import matplotlib.pyplot as plt
from CAPTrajectory import CAPHamiltonian

#Change these lines to suit your system
##########################################
ref_energy = -109.36009153
guess = 2.2
eta_list = np.linspace(0,500,101)
eta_list = eta_list * 1E-5
RASSI_FILE = "../../opencap/openmolcas/nosymm.rassi.h5"
OUTPUT_FILE = "../../opencap/openmolcas/nosymm.out"
nbasis = 119
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

f = h5py.File(RASSI_FILE, 'r')
arr = f["SFS_TRANSITION_DENSITIES"]
arr2 = np.array(f["AO_OVERLAP_MATRIX"])
arr2 = np.reshape(arr2,(nbasis,nbasis))

# Method 1: Read data from output
s = pyopencap.System(sys_dict)
s.check_overlap_mat(arr2,"openmolcas",RASSI_FILE)
pc = pyopencap.CAP(s,cap_dict,10,"openmolcas")
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()
h0 = pc.get_H()


# Method 2: Separate alpha and beta densities
pc = pyopencap.CAP(s,cap_dict,10,"openmolcas")
for i in range(0,10):
    for j in range(i,10):
        arr1 = 0.5*np.reshape(arr[i][j],(nbasis,nbasis))
        pc.add_tdms(arr1,arr1,i,j,"openmolcas",RASSI_FILE)
        if i!=j:
            pc.add_tdms(arr1,arr1,j,i,"openmolcas",RASSI_FILE)
pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()


# Method 3: Spin traced densities
pc = pyopencap.CAP(s,cap_dict,10,"openmolcas")
for i in range(0,10):
    for j in range(i,10):
        arr1 = np.reshape(arr[i][j],(nbasis,nbasis))
        pc.add_tdm(arr1,i,j,"openmolcas",RASSI_FILE)
        if i!=j:
            pc.add_tdm(arr1,j,i,"openmolcas",RASSI_FILE)
pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()
print(DataFrame(h0).to_string(index=False, header=False))

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

