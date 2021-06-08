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
CAP/FCI on doubly excited H2 Feshbach resonance. Basis set is:
'''

import pyopencap
from pandas import DataFrame
from pyscf import gto, scf, ci, ao2mo, tools, fci
import numpy as np
import time

start = time.time()


cap_dict = {
            "cap_type": "box",
            "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
            "Radial_precision": "14",
            "angular_points": "110"
}
nstates = 20

# pyscf mol
ghost_bas = gto.basis.load('h2.nw', 'N')
H_bas = gto.basis.load('h2.nw', 'H')
mol = gto.M(
            atom = 'ghost 0.0000000000 0.0000000000 -0.7400000000; \
            ghost 0.0000000000 0.0000000000 0.7400000000;          \
            H 0.0000000000 0.0000000000 0.3705000000;           \
             H 0.0000000000 0.0000000000 -0.3705000000',
             basis = {'H': H_bas, 'ghost': ghost_bas},
             symmetry='d2h')
mol.build()
myhf = scf.RHF(mol)
myhf.kernel()

# create system using molden file
molden_dict = {"basis_file":"molden_in.molden",
    "molecule": "molden"}
tools.molden.from_scf(myhf,"molden_in.molden")
s = pyopencap.System(molden_dict)
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")

# run full ci
fs = fci.FCI(mol, myhf.mo_coeff)
fs.nroots = nstates
e, c = fs.kernel()

# fill density matrices
h0 = np.zeros((nstates,nstates))
for i in range(0,len(e)):
    h0[i][i]=e[i]

# compute CAP
pc = pyopencap.CAP(s,cap_dict,nstates)
pc.compute_ao_cap()
for i in range(0,len(fs.ci)):
    for j in range(0,len(fs.ci)):
        dm1 = fs.trans_rdm1(fs.ci[i],fs.ci[j],myhf.mo_coeff.shape[1],mol.nelec)
        dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
        pc.add_tdm(dm1_ao,i,j,"pyscf")


pc.compute_projected_cap()
W=pc.get_projected_cap()

from pyopencap.analysis import CAPHamiltonian as CAPH
my_CAPH = CAPH(H0=h0,W=W)

my_CAPH.export("H2_fci_opencap.out")

end = time.time()
print("Time:")
print(end - start)

import os
os.remove("molden_in.molden")
