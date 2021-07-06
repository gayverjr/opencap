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
from adcc.State2States import State2States
import pyopencap
import psi4
from pandas import DataFrame
import os
import time

start = time.time()
psi4.core.set_output_file('output.dat', False)

mol = psi4.geometry("""
    N 0.0000000000 0.0000000000 0.548756750
    N 0.0000000000 0.0000000000 -0.548756750
    Gh(He) 0.0000000000 0.0000000000 0.000
    Symmetry C1
    """)


def basisspec_psi4_yo__anonymous03952cbd(mol, role):
    basstrings = {}
    mol.set_basis_all_atoms("test_bas", role=role)
    basstrings['test_bas'] = """
        spherical
        ****
        N     0
        S    10   1.00
              1.142000D+04           5.230000D-04
              1.712000D+03           4.045000D-03
              3.893000D+02           2.077500D-02
              1.100000D+02           8.072700D-02
              3.557000D+01           2.330740D-01
              1.254000D+01           4.335010D-01
              4.644000D+00           3.474720D-01
              1.293000D+00           4.126200D-02
              5.118000D-01          -8.508000D-03
              1.787000D-01           2.384000D-03
        S    10   1.00
              1.142000D+04          -1.150000D-04
              1.712000D+03          -8.950000D-04
              3.893000D+02          -4.624000D-03
              1.100000D+02          -1.852800D-02
              3.557000D+01          -5.733900D-02
              1.254000D+01          -1.320760D-01
              4.644000D+00          -1.725100D-01
              1.293000D+00           1.518140D-01
              5.118000D-01           5.999440D-01
              1.787000D-01           3.874620D-01
        S    1   1.00
              1.293000D+00           1.000000D+00
        S    1   1.00
              1.787000D-01           1.000000D+00
        S    1   1.00
              0.0576000              1.0000000
        P    5   1.00
              2.663000D+01           1.467000D-02
              5.948000D+00           9.176400D-02
              1.742000D+00           2.986830D-01
              5.550000D-01           4.984870D-01
              1.725000D-01           3.370230D-01
        P    1   1.00
              5.550000D-01           1.000000D+00
        P    1   1.00
              1.725000D-01           1.000000D+00
        P    1   1.00
              0.0491000              1.0000000
        D    1   1.00
              1.654000D+00           1.000000D+00
        D    1   1.00
              4.690000D-01           1.000000D+00
        D    1   1.00
              0.1510000              1.0000000
        F    1   1.00
              1.093000D+00           1.0000000
        F    1   1.00
              0.3640000              1.0000000
        ****
        He     0
        S   1   1.00
        0.0288000              1.0000000
        S   1   1.00
        0.0144000              1.0000000
        S   1   1.00
        0.0072000              1.0000000
        S   1   1.00
        1.00E-08               1.0000000
        P   1   1.00
        0.0245500              1.0000000
        P   1   1.00
        0.0122750              1.0000000
        P   1   1.00
        0.0061375              1.0000000
        D   1   1.00
        0.0755000              1.0000000
        D   1   1.00
        0.0377500              1.0000000
        D   1   1.00
        0.0188750              1.0000000
        ****
        """
    return basstrings

psi4.qcdb.libmintsbasisset.basishorde['ANONYMOUS03952CBD'] = basisspec_psi4_yo__anonymous03952cbd
psi4.core.set_global_option("BASIS", "anonymous03952cbd")
E, wfn = psi4.energy('scf', return_wfn=True)
psi4.molden(wfn, 'n2.molden')
# add 7F to molden file, psi4 doesn't write it for some reason
with open("n2.molden", "a") as myfile:
    myfile.write("\n [7F] \n")


# checking overlap matrix
mints = psi4.core.MintsHelper(wfn.basisset())
S_mat = np.asarray(mints.ao_overlap())
molden_dict = {"basis_file":"n2.molden",
    "molecule": "molden"}
s = pyopencap.System(molden_dict)
s.check_overlap_mat(S_mat,"psi4")
print('', end='',flush=True)
cap_dict = {
    "cap_type": "box",
        "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
}
nstates = 30
pc = pyopencap.CAP(s,cap_dict,nstates)

# ground state energy
state = adcc.adc2(wfn, n_singlets=1)
E_0 = E + state.ground_state.energy_correction(2)
print("MP2 energy of neutral:" + str(E_0))

# now starting anion calculation, starting from UHF reference
mol = psi4.geometry("""
    -1 2
    N 0.0000000000 0.0000000000 0.548756750
    N 0.0000000000 0.0000000000 -0.548756750
    Gh(He) 0.0000000000 0.0000000000 0.000
    Symmetry C1
    """)
psi4.qcdb.libmintsbasisset.basishorde['ANONYMOUS03952CBD'] = basisspec_psi4_yo__anonymous03952cbd
psi4.core.set_global_option("BASIS", "anonymous03952cbd")
psi4.core.set_global_option("REFERENCE", "uhf")
E, wfn = psi4.energy('scf', return_wfn=True)

state = adcc.adc2(wfn, n_states=nstates)
print(state.describe())
E_0 = E + state.ground_state.energy_correction(2)
print("MP2 energy of anion:" + str(E_0))

h0 = np.zeros((nstates,nstates))
for i in range(0,len(h0)):
    h0[i][i]= E_0 + state.excitation_energy[i]

for i in range(0,nstates):
    state2state = State2States(state, initial=i)
    rdm = state.state_dm[i]
    rdm_alpha,rdm_beta = rdm.to_ao_basis()
    rdm_alpha = rdm_alpha.to_ndarray()
    rdm_beta = rdm_beta.to_ndarray()
    pc.add_tdms(rdm_alpha,rdm_beta,i,i,"psi4")
    for j,tdm in enumerate(state2state.transition_dm):
        tdm_alpha,tdm_beta = tdm.to_ao_basis()
        tdm_alpha = tdm_alpha.to_ndarray()
        tdm_beta = tdm_beta.to_ndarray()
        pc.add_tdms(tdm_alpha,tdm_beta,i,i+j+1,"psi4")
        pc.add_tdms(tdm_alpha.conj().T,tdm_beta.conj().T,i+j+1,i,"psi4")

pc.compute_projected_cap()
W = pc.get_projected_cap()

from pyopencap.analysis import CAPHamiltonian as CAPH
my_CAPH = CAPH(H0=h0,W=W)

my_CAPH.export("n2_adc_opencap.out")

end = time.time()
print("Time:")
print(end - start)
