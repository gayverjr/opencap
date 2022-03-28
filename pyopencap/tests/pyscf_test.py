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
    SOFTWARE.'''

import pyopencap
import numpy as np
import os
import sys
import pytest
try:
    import pyscf
    from pyscf import gto, scf, fci, tools
except:
    pass

sys_dict = {"molecule": "molden", "basis_file": "molden_file.molden"}

cap_dict = {
    "cap_type": "box",
    "cap_x": "6.00",
    "cap_y": "6.00",
    "cap_z": "6.7",
    "Radial_precision": "14",
    "angular_points": "110"
}


def get_hf():
    mol = gto.M(atom='H        0.0     0.0     0.54981512; \
                H        0.0     0.0     -0.54981512;',
                basis="cc-pvdz")
    mol.build()
    myhf = scf.RHF(mol)
    myhf.kernel()
    return mol, myhf


@pytest.mark.skipif('pyscf' not in sys.modules,
                    reason="requires the pyscf library")
def test_from_molden():
    mol, myhf = get_hf()
    tools.molden.from_scf(myhf, "molden_file.molden")
    s = pyopencap.System(sys_dict)
    pyscf_smat = scf.hf.get_ovlp(mol)
    s.check_overlap_mat(pyscf_smat, "pyscf")
    os.remove("molden_file.molden")

@pytest.mark.skipif('pyscf' not in sys.modules,
                    reason="requires the pyscf library")
def test_ovlp_ordering():
    mol, myhf = get_hf()
    tools.molden.from_scf(myhf, "molden_file.molden")
    s = pyopencap.System(sys_dict)
    pyscf_smat = scf.hf.get_ovlp(mol)
    os.remove("molden_file.molden")
    assert np.allclose(s.get_overlap_mat("pyscf"),pyscf_smat)


@pytest.mark.skipif('pyscf' not in sys.modules,
                    reason="requires the pyscf library")
def test_pyscf():
    mol, myhf = get_hf()
    tools.molden.from_scf(myhf, "molden_file.molden")
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s, cap_dict, 3)
    pc.compute_ao_cap(cap_dict)
    W_AO = pc.get_ao_cap("pyscf")
    W2 = np.zeros((3,3))
    fs = fci.FCI(mol, myhf.mo_coeff)
    fs.nroots = 3
    e, c = fs.kernel()
    for i in range(0, len(fs.ci)):
        for j in range(0, len(fs.ci)):
            dm1 = fs.trans_rdm1(fs.ci[i], fs.ci[j], myhf.mo_coeff.shape[1],
                                mol.nelec)
            dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1,
                               myhf.mo_coeff.conj())
            pc.add_tdm(dm1_ao, i, j, "pyscf")
            W2[i,j] = -1.0 * np.trace(W_AO,dm1_ao)
    pc.compute_projected_cap()
    W = pc.get_projected_cap()
    os.remove("molden_file.molden")
    assert np.allclose(W2,W)
