'''Copyright (c) 2021 James Gayvert
    
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

'''
CAP/ADC(2) on N2-. Basis set is aug-cc-pvtz + 3s3p3d[Gh]. Requires ADCC.
See https://adc-connect.org/
'''

import adcc
import numpy as np
from pyscf import gto, scf,tools
from adcc.State2States import State2States
import pyopencap
import time

start = time.time()

sys_dict = {"geometry":'''    N 0.0000000000 0.0000000000 0.548756750
    N 0.0000000000 0.0000000000 -0.548756750
    X       0.0     0.0     0.0''',
        "basis_file":"n2.bas",
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
pc = pyopencap.CAP(s,cap_dict,nstates)


ghost_bas = gto.basis.load('n2.nw', 'X')
N_bas = gto.basis.load('n2.nw', 'N')
mol = gto.M(
            atom = '    N 0.0000000000 0.0000000000 0.548756750;\
            N 0.0000000000 0.0000000000 -0.548756750;\
            ghost      0.0     0.0     0.0',
            basis = {'N': N_bas, 'ghost': ghost_bas} )
mol.build()
scfres = scf.RHF(mol)
scfres.kernel()

# check overlap matrix
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")

# run scf and mp2 on neutral
state = adcc.adc2(scfres, n_singlets=1)
E_0 = scfres.energy_tot() + state.ground_state.energy_correction(2)
print("MP2 energy of neutral:" + str(E_0))

# run scf and mp2 on anion
mol.charge = -1
mol.spin = 1
mol.build()
scfres = scf.UHF(mol)
scfres.kernel()

# now run adc
state = adcc.adc2(scfres, n_states=nstates)
print(state.describe())
E_0 = scfres.energy_tot() + state.ground_state.energy_correction(2)
print("MP2 energy of anion:" + str(E_0))

# define H0
h0 = np.zeros((nstates,nstates))
for i in range(0,len(h0)):
    h0[i][i]= E_0 + state.excitation_energy[i]

# calculate W
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

from pyopencap.analysis import CAPHamiltonian as CAPH
my_CAPH = CAPH(H0=h0,W=W)

my_CAPH.export("n2_adc_opencap.out")

end = time.time()
print("Time:")
print(end - start)

