import numpy as np
import matplotlib.pyplot as plt
from CAPTrajectory import CAPHamiltonian

#########################################
# Alter these values to suit your purposes
ref_energy = 0
guess = 2.5
eta_list = np.linspace(0,5000,101)
###########################
au2eV= 27.2113961
eta_list = eta_list * 1E-5

CAPH = CAPHamiltonian(output_file = "n2_opencap.out")

CAPH.run_trajectory(eta_list)

traj = CAPH.track_state(1,tracking="overlap")
uc_energy,uc_eta_opt = traj.find_eta_opt()
corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

print((uc_energy-ref_energy)*27.2114)
print((corr_energy-ref_energy)*27.2114)

plt.plot(np.real(traj.uncorrected_energies),np.imag(traj.uncorrected_energies),'-ro',label='Uncorrected Trajectory')
plt.plot(np.real(traj.corrected_energies),np.imag(traj.corrected_energies),'-bo',label='Corrected Trajectory')
plt.legend()
plt.show()




