import adcc
import numpy as np
from pyscf import gto, scf,tools
from adcc.State2States import State2States
import pyopencap
from pandas import DataFrame
import time

start = time.time()

sys_dict = {"geometry":'''N        0.0     0.0     0.54981512
    N        0.0     0.0     -0.54981512
    X       0.0     0.0     0.0''',
        "basis_file":"n2_psi4format.bas",
        "molecule": "inline"
}

s = pyopencap.System(sys_dict)

cap_dict = {
    "cap_type": "box",
    "cap_x":"2.76",
    "cap_y":"2.76",
    "cap_z":"4.88",
}

nstates = 30
pc = pyopencap.CAP(s,cap_dict,nstates,"pyscf")


ghost_bas = gto.basis.load('n2_nwchemformat.bas', 'X')
N_bas = gto.basis.load('n2_nwchemformat.bas', 'N')
mol = gto.M(
            atom = 'N        0.0     0.0     0.54981512; \
            N        0.0     0.0     -0.54981512;        \
            ghost      0.0     0.0     0.0',
            basis = {'N': N_bas, 'ghost': ghost_bas} )
mol.build()
scfres = scf.RHF(mol)
scfres.kernel()
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")
state = adcc.adc2(scfres, n_singlets=1)
E_0 = scfres.energy_tot() + state.ground_state.energy_correction(2)
print("E_0:" + str(E_0))

mol.charge = -1
mol.spin = 1
mol.build()
scfres = scf.UHF(mol)
scfres.kernel()
state = adcc.adc2(scfres, n_states=nstates)
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


pc.compute_projected_cap()
W=pc.get_projected_cap()

print("Printing out matrices required for Projected CAP calculation.")
print("Number of states: " + str(nstates))
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP Matrix")
print(DataFrame(W).to_string(index=False, header=False))

end = time.time()
print("Time:")
print(end - start)

