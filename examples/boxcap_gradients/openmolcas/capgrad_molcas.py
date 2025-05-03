#!/usr/bin/env python
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
Evaluating cuboid CAP gradients with pre-computed energies and 1-TDMs!

XMS-CASPT2 on dinitrogen anion (with ghost atom in COM) with aug-cc-pvtz+X[3s3p3d] basis set.
'''


import os, sys
import scipy.linalg as LA
from functools import reduce
import numpy as np
import pandas as pd

opencap_driver = False
try:
    __import__("pyopencap")
    print("\nOpenCAP found and will be called in %s \n" % (sys.argv[0]))
    opencap_driver = True
except ImportError:
    print("Not found pyopencap")


# Define some variables
allroots = 5 #Total roots in state-averaging
eta_opt = 0.009 # Fixed eta value
track_root = 2 # 2nd root (starting 1) is the resonance
natom = 3 # Including one ghost atom!


# Some utility function
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

def parse_PCAP_HAMILTONIAN(RASSI_FILE, OUTPUT_FILE):
    import pyopencap
    from pyopencap.analysis import CAPHamiltonian


    sys_dict = {"molecule": "molcas_rassi", "basis_file": RASSI_FILE}
    s = pyopencap.System(sys_dict)
    cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88", "do_numerical" : "True",
            }

    es_dict = { "package": "openmolcas",
    "method" : "ms-caspt2",
        "molcas_output":OUTPUT_FILE ,
            "rassi_h5":RASSI_FILE}

    pc = pyopencap.CAP(s, cap_dict, allroots)


    pc.read_data(es_dict)
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    H0 = pc.get_H()
    CAPH = CAPHamiltonian(H0=H0,W=W)

    CAPH.export("./CAPMAT.out")

    return pc, CAPH, H0, W


# Generate projected CAP object, Hamiltonians with openmolcas generated TDM and energy files.
OUTPUT_FILE = "../../analysis/N2/ref_outputs/xms.out"
RASSI_FILE = "../../analysis/N2/ref_outputs/xms.rassi.h5"
pc, CAPH, H0, W = parse_PCAP_HAMILTONIAN(RASSI_FILE, OUTPUT_FILE)

# Generate rotation vectors for a fixed eta
eigv, Leigvc, Reigvc = run_traj(H0, W, eta_opt)


# Call all derivative terms (cuboid CAP only)
pc.compute_projected_capG()
_WG = pc.get_projected_capG()
pc.compute_projected_cap_der()
_WD = pc.get_projected_cap_der()
# All opencap calculations done

# Print them in a nice (some say unnecessary!) fashion
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

for i, state1 in enumerate([track_root]):
    state1 -= 1
    # Storing only real part of the gradients
    QMout['grad_diag_res'].append([[grad_store[_atom]['x'][state1].real,
                                    grad_store[_atom]['y'][state1].real,
                                    grad_store[_atom]['z'][state1].real] for _atom in range(natom)])


float_format = lambda x: f"{x: 16.12f}"

for _, (state_gradients, state) in enumerate(zip(QMout['grad_diag_res'], [track_root])):
    print(f"\nState: {state} (hartree/bohr)")
    df = pd.DataFrame(state_gradients, columns=['X', 'Y', 'Z'])
    df.insert(0, "Atom no", range(1, len(df) + 1))
    df.iloc[:, 1:] = df.iloc[:, 1:].applymap(float_format)
    print(df.to_string(index=False))
