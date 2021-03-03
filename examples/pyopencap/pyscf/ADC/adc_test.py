import pyopencap
import numpy as np
from pandas import DataFrame
import h5py
import matplotlib.pyplot as plt
from CAPTrajectory import CAPHamiltonian

ref_energy = -109.39504305049073
eta_list = np.linspace(0,5000,101)
eta_list*=1E-5

CAPH = CAPHamiltonian(output_file="adc_example.out")

CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.energies_ev(ref_energy=ref_energy)),np.imag(CAPH.energies_ev(ref_energy=ref_energy)),'ro',label='Uncorrected Trajectory')
plt.show()

traj = CAPH.track_state(12,tracking="overlap")
uc_energy,uc_eta_opt = traj.find_eta_opt()
corr_energy,corr_eta_opt = traj.find_eta_opt(corrected=True)

print((uc_energy-ref_energy)*27.2114)
print((corr_energy-ref_energy)*27.2114)

plt.plot(np.real(traj.energies_ev(ref_energy=ref_energy)),np.imag(traj.energies_ev(ref_energy=ref_energy)),'-ro',label='Uncorrected Trajectory')
plt.plot(np.real(traj.energies_ev(ref_energy=ref_energy,corrected=True)),np.imag(traj.energies_ev(ref_energy=ref_energy,corrected=True)),'-bo',label='Corrected Trajectory')
plt.legend()
plt.show()
