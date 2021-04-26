import adcc
import numpy as np
from pyscf import gto, scf,tools
from adcc.State2States import State2States
import pyopencap
from pandas import DataFrame
import matplotlib.pyplot as plt
from CAPTrajectory import CAPHamiltonian
molden_dict = {"basis_file":"basis.molden",
    "molecule": "molden"}
s = pyopencap.System(molden_dict)

cap_dict = {
    "cap_type": "box",
        "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

nstates = 30
pc = pyopencap.CAP(s,cap_dict,nstates,"pyscf")


mol=tools.molden.load('basis.molden')[0]
mol.build()
scfres = scf.RHF(mol)
scfres.kernel()
state = adcc.adc3(scfres, n_singlets=1)
E_0 = scfres.energy_tot() + state.ground_state.energy_correction(2)
print("E_0:" + str(E_0))


mol=tools.molden.load('basis.molden')[0]
mol.charge = -1
mol.spin = 1
mol.build()
scfres = scf.UHF(mol)
scfres.kernel()
state = adcc.adc2(scfres, n_states=nstates)
#E_0 = scfres.energy_tot() + state.ground_state.energy_correction(2)
print(state.describe())

h0 = np.zeros((nstates,nstates))
for i in range(0,len(h0)):
    h0[i][i]= E_0 + state.excitation_energy[i]

for i in range(0,nstates):
    state2state = State2States(state, initial=i)
    rdm = state.state_dm[i]
    rdm_alpha,rdm_beta = rdm.to_ao_basis()
    rdm_alpha = rdm_alpha.to_ndarray()
    rdm_beta = rdm_beta.to_ndarray()
    pc.add_tdms(rdm_alpha,rdm_beta,i,i,"pyscf")
    for j,tdm in enumerate(state2state.transition_dm):
        tdm_alpha,tdm_beta = tdm.to_ao_basis()
        tdm_alpha = tdm_alpha.to_ndarray()
        tdm_beta = tdm_beta.to_ndarray()
        pc.add_tdms(tdm_alpha,tdm_beta,i,i+j+1,"pyscf")
        pc.add_tdms(tdm_alpha,tdm_beta,i+j+1,i,"pyscf")

pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()

("Printing out matrices required for Perturbative CAP calculation.")
print("Number of states: " + str(nstates))
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP Matrix")
print(DataFrame(mat).to_string(index=False, header=False))

eta_list = np.linspace(0,5000,101)

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
