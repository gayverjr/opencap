import adcc
import numpy as np
from pyscf import gto, scf,tools
from adcc.State2States import State2States
import pyopencap
from pandas import DataFrame
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

nstates = 50
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
    print("Added state DM for state:" + str(i))
    for j,tdm in enumerate(state2state.transition_dm):
        tdm_alpha,tdm_beta = tdm.to_ao_basis()
        tdm_alpha = tdm_alpha.to_ndarray()
        tdm_beta = tdm_beta.to_ndarray()
        print("Added transition dms for states:" + str(i) + "," + str(i+j+1))
        pc.add_tdms(tdm_alpha,tdm_beta,i,i+j+1,"pyscf")
        pc.add_tdms(tdm_alpha,tdm_beta,i+j+1,i,"pyscf")


print("Calculating cap matrix")
pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()
("Printing out matrices reiquired for Perturbative CAP calculation.")
print("Number of states: " + str(nstates))
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP Matrix")
print(DataFrame(mat).to_string(index=False, header=False))
