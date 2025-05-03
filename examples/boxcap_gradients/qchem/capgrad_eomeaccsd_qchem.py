'''Copyright (c) 2024 Soubhik Mondal

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
Evaluating cuboid CAP gradients with pre-computed energies and 1-TDMs.

pCAP-EOM-EA-CCSD (5 EA states) on N2- with aug-cc-pvtz+X[3s3p3d] basis set.

 State         Re(E)            Im(E)
   1      2.629737516514  -0.193944247865

'''

import numpy as np
import scipy.linalg as LA
from functools import reduce
import pyopencap
import pandas as pd
EV_TO_AU = 1/27.21138602

# Some utility functions
def _sort_eigenvectors(eigv, eigvc):
    '''
    Sort eigen vector columns according to the ascending order of
    real part of eigen values.
    :param eigv:
    :param eigvc:
    :return: eigv, eigvc
    '''
    idx = eigv.real.argsort()
    eigv = eigv[idx]
    eigvc = eigvc[:, idx]
    return eigv, eigvc

def run_traj(H0, W, eta_opt):
    CAPHAM = H0 + 1.0j * eta_opt * W
    CAPHAM = 0.5 * (CAPHAM.T + CAPHAM)
    eigv, Reigvc = _sort_eigenvectors(*LA.eig(CAPHAM))
    W = reduce(np.dot, (Reigvc.T, Reigvc))
    W_sqrt = LA.sqrtm(W)
    W_inv_sqrt = LA.inv(W_sqrt)
    Reigvc = reduce(np.dot, (Reigvc, W_inv_sqrt))
    return eigv, Reigvc, Reigvc

def get_ref_energy(finame):
    with open(finame,'r') as f:
        lines = f.readlines()
    for l in lines:
        if 'CCSD total energy' in l:
            return float(l.split()[-1])




nstates = 5 # Total 5 EA states
eta_opt = 0.0075
natom = 3
track_root = [1]

# Files to parse information from
FCHK_FILE = "../../analysis/N2/ref_outputs/qc_inp.fchk"
OUTPUT_FILE = "../../analysis/N2/ref_outputs/qc_inp.out"


# Get CCSD energy
ref_energy = get_ref_energy(OUTPUT_FILE)
print('Ref energy: ', ref_energy)

# Build opencap dictionaries
sys_dict = {"molecule": "qchem_fchk", "basis_file": FCHK_FILE}
es_dict = {
    "package": "qchem",
    "method": "eom",
    "qchem_output": OUTPUT_FILE,
    "qchem_fchk": FCHK_FILE,
}

cap_dict = {"cap_type": "box",
            "cap_x": "2.760",
            "cap_y": "2.760",
            "cap_z": "4.880",
            "Radial_precision": "16",
            "angular_points": "590",
            "do_numerical": "true"
}

# Get pyopencap object
sys = pyopencap.System(sys_dict)
pc = pyopencap.CAP(sys, cap_dict, nstates)
pc.read_data(es_dict)
pc.compute_projected_cap()
W = pc.get_projected_cap()
H0 = pc.get_H()


# Get rotation vectors for a fixed eta

eigv, Leigvc, Reigvc = run_traj(H0, W, eta_opt)


float_format = lambda x: f"{x: 16.12f}"
print(f"\n\nState energies:(eV)")

df = pd.DataFrame({
    "EA root": range(1, len(eigv) + 1),
    "Re(E)": [(ev.real - ref_energy) / EV_TO_AU for ev in eigv],
    "Im(E)": [ev.imag / EV_TO_AU for ev in eigv]
})

df.iloc[:, 1:] = df.iloc[:, 1:].applymap(float_format)

print(df.to_string(index=False))
print('\n')


# Call all derivative terms (cuboid CAP only)
pc.compute_projected_capG()
_WG = pc.get_projected_capG()
pc.compute_projected_cap_der()
_WD = pc.get_projected_cap_der()
# All opencap calculations done



# Print them in a nice fashion
Grad_correct = []
for iatom in range(natom):
    Grad_correct.append({'x': _WG[iatom]['x']+_WD[iatom]['x'],
                         'y': _WG[iatom]['y']+_WD[iatom]['y'],
                         'z': _WG[iatom]['z']+_WD[iatom]['z']})

QMout = {}
QMout['grad_diag_res'] = []
grad_store = []
for j in range(natom):
    grad_diag = {}
    for cart in ['x', 'y', 'z']:
        _grad_MCH = 1.0j * eta_opt * Grad_correct[j][cart]
        _diag = Leigvc.T @ _grad_MCH @ Reigvc
        grad_diag[cart] = np.diagonal(_diag).real

    grad_store.append({'x': grad_diag['x'],
                       'y': grad_diag['y'],
                       'z': grad_diag['z']})

for i, state1 in enumerate(track_root):
    state1 -= 1
    # Storing only real part of the gradients
    QMout['grad_diag_res'].append([[grad_store[_atom]['x'][state1].real,
                                    grad_store[_atom]['y'][state1].real,
                                    grad_store[_atom]['z'][state1].real] for _atom in range(natom)])


float_format = lambda x: f"{x: 16.12f}"

for _, (state_gradients, state) in enumerate(zip(QMout['grad_diag_res'], track_root)):
    print(f"\nState: {state} (hartree/bohr)")
    df = pd.DataFrame(state_gradients, columns=['X', 'Y', 'Z'])
    df.insert(0, "Atom no", range(1, len(df) + 1))
    df.iloc[:, 1:] = df.iloc[:, 1:].applymap(float_format)
    print(df.to_string(index=False))
