#!/usr/bin/env python

import os, sys
import scipy.linalg as LA
from functools import reduce
import numpy as np
from tabulate import tabulate

opencap_driver = False
try:
    __import__("pyopencap")
    print("\nOpenCAP found and will be called in %s \n" % (sys.argv[0]))
    opencap_driver = True
except ImportError:
    print("Not found pyopencap")


# Define some variables
EV_TO_HARTREE = 1/27.211386245988
allroots=23 #Total roots in state-averaging
fname="./ref_outputs/iodoethene/MOLCAS"
eta_opt = 2128.0E-5 # Fixed eta value
track_root = 10 # 10th root (starting 1) is the resonance
natom = 6 # Including one ghost atom!


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
    if opencap_driver:
        os.environ['OMP_NUM_THREADS']=str(10)

    import pyopencap
    from pyopencap.analysis import CAPHamiltonian


    sys_dict = {"molecule": "molcas_rassi", "basis_file": RASSI_FILE}
    s = pyopencap.System(sys_dict)
    cap_dict = {"cap_type": "box",
                "cap_x":"18.55",
                "cap_y":"7.36",
                "cap_z":"11.70","Radial_precision": "16",
                "angular_points": "590", "do_numerical": "true"}

    pc = pyopencap.CAP(s, cap_dict, allroots)

    es_dict = {"package": "openmolcas",
                "method": "sa-casscf",
                "molcas_output": OUTPUT_FILE,
                "rassi_h5": RASSI_FILE}

    pc.read_data(es_dict)
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    H0 = pc.get_H()
    CAPH = CAPHamiltonian(H0=H0,W=W)

    CAPH.export("./CAPMAT.out")

    return pc, CAPH, H0, W


# Generate projected CAP object, Hamiltonians with openmolcas generated TDM and energy files.
RASSI_FILE = "%s.rassi.h5" % fname
OUTPUT_FILE = "%s.out" % fname
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


float_format = lambda x: "{: 14.8f}".format(x)
for state_index, state_gradients in enumerate(QMout['grad_diag_res'], start=1):
    print(f"State: {state_index}")
    headers = ['Atom no', 'X', 'Y', 'Z']
    table_data = []
    for atom_index, atom_gradients in enumerate(state_gradients, start=1):
        formatted_gradients = [float_format(grad) for grad in atom_gradients]
        table_data.append([atom_index] + formatted_gradients)
    print(tabulate(table_data, headers=headers, tablefmt="pretty", stralign="right"))

