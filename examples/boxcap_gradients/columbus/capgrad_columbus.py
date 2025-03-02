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

Example 1: MR-CISD on N2- with aug-cc-pvtz+3s3p3d basis set.
Example 2: SA-MCSCF on N2- with cc-pvtz+2s5p2d basis set.
'''



import pyopencap
import numpy as np
import scipy.linalg as LA
from pyopencap.analysis import colparser, colparser_mc
from functools import reduce
import pandas as pd


# Some utility functions
def _sort_eigenvectors(eigv, eigvc):
    '''
    Sort eigen vector columns according to the ascending order of
    real part of eigen values.
    :param np.ndarray eigv: eigenvalues
    :param  np.ndarray eigvc: eigen vectors
    :return: sorted eigen values and vectors
    '''
    idx = eigv.real.argsort()
    eigv = eigv[idx]
    eigvc = eigvc[:, idx]
    return eigv, eigvc

def run_traj(H0, W, eta_opt):
    '''
    For a given CAP augmented H^CAP matrix, it performs diagonalization!
    :param H0: Zero Order energy matrix
    :param W: CAP matrix
    :param eta_opt: A fixed eta value
    :return: eigv, Reigvc, Leigvc: complex eigen-values, rotation vectors (left, right)
    '''
    CAPHAM = H0 + 1.0j * eta_opt * W
    CAPHAM = 0.5 * (CAPHAM.T + CAPHAM)
    eigv, Reigvc = _sort_eigenvectors(*LA.eig(CAPHAM))
    W = reduce(np.dot, (Reigvc.T, Reigvc))
    W_sqrt = LA.sqrtm(W)
    W_inv_sqrt = LA.inv(W_sqrt)
    Reigvc = reduce(np.dot, (Reigvc, W_inv_sqrt))
    return eigv, Reigvc, Reigvc




'''
Example 1: MR-CISD on N2- with aug-cc-pvtz+3s3p3d basis set.
'''


print("\n\n*** MR-CISD on N2- with aug-cc-pvtz+3s3p3d basis set.\n\n")
# Define some variables for N2 molecule
nstates = 4 # SA-4 calculations
natom = 3 #Two atoms of N + 1 ghost atom at COM
eta_opt = 0.00848
track_root = [1, 2, 3, 4]
# CAP gradients is calculated for all states by-default
# In this example: 2nd state is resonance state


# Define CAP-INPUTS: energies and TDMs
path_TO_CAPINP = "../../pyopencap/columbus/data_files/"
parser = colparser('%s/molden_mo_mc.sp'%path_TO_CAPINP, f"{path_TO_CAPINP}/tranls")
molden_dict = {"basis_file":"%s/molden_mo_mc.sp"%path_TO_CAPINP, "molecule": "molden"}
s = pyopencap.System(molden_dict)
#

# Set up CAP calculations with box-CAP
cap_dict = {"cap_type": "box",
            "cap_x":"2.77",
            "cap_y":"2.77",
            "cap_z":"4.88",
            "Radial_precision": "16",
            "angular_points": "590",
            "do_numerical": "true"}

pc = pyopencap.CAP(s,cap_dict,nstates)

# Build W matrix
for i in range(0,nstates):
    for j in range(i,nstates):
        if i==j:
                dm1_ao = parser.sdm_ao(i+1, filename=f"{path_TO_CAPINP}/cid1trfl.FROMdrt1.state{i+1}TOdrt1.state{i+1}.iwfmt")
                pc.add_tdm(dm1_ao, i, j, 'molden')
        else:
                dm1_ao = parser.tdm_ao(i+1, j+1, filename=f"{path_TO_CAPINP}/cid1trfl.FROMdrt1.state{i+1}TOdrt1.state{j+1}.iwfmt")
                pc.add_tdm(dm1_ao,i,j,'molden')
                pc.add_tdm(dm1_ao.T.conj(),j,i,'molden')
pc.compute_projected_cap()
W=pc.get_projected_cap()

# Generate rotation vectors for a fixed eta
H0 = parser.get_H0(filename='%s/ciudgsm'%path_TO_CAPINP, correction_type='eci')
# Only available for energy type eci, other corrections are apriori
eigv, Leigvc, Reigvc = run_traj(H0, W, eta_opt)


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


'''
Example 2: SA-MCSCF on N2- with cc-pvtz+2s5p2d basis set.
'''

print("\n\n*** SA-MCSCF on N2- with cc-pvtz+2s5p2d basis set.\n\n")

# Define some variables for N2 molecule
nstates = 9 # SA-9 calculations
natom = 2 #Two atoms of N
eta_opt = 0.0079
track_root = [2]
# In this example: 2nd state is resonance state


# Define CAP-INPUTS: energies and TDMs (MCSCF)
path_TO_CAPINP = "../../pyopencap/columbus/data_files_mcscf/"
parser = colparser_mc('%s/molden_mo_mc.sp'%path_TO_CAPINP)
molden_dict = {"basis_file":"%s/molden_mo_mc.sp"%path_TO_CAPINP, "molecule": "molden"}
s = pyopencap.System(molden_dict)
# CAP_dict is same as last one
pc = pyopencap.CAP(s,cap_dict,nstates)

# Build W matrix in MCSCF basis
for i in range(0,nstates):
    for j in range(i,nstates):
        if i==j:
                fnameIN = "%s/mcsd1fl.drt1.st%2.2i.iwfmt"%(path_TO_CAPINP, i+1)
                dm1_ao = parser.sdm_ao(i+1,filename=fnameIN)
                pc.add_tdm(dm1_ao,i,j,'molden')
        else:
                fnameIN ="%s/mcsd1fl.drt1.st%2.2i-st%2.2i.iwfmt"%(path_TO_CAPINP, i+1, j+1)
                dm1_ao = parser.tdm_ao(i+1, j+1, filename=fnameIN)
                pc.add_tdm(dm1_ao,i,j,'molden')
                pc.add_tdm(dm1_ao.T,j,i,'molden')
pc.compute_projected_cap()
W=pc.get_projected_cap()

# Generate rotation vectors for a fixed eta
H0 = parser.get_H0(filename='%s/mcscfsm'%path_TO_CAPINP,)
eigv, Leigvc, Reigvc = run_traj(H0, W, eta_opt)


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
