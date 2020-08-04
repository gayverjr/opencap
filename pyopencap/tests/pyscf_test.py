import pyopencap as pycap
from pyscf import gto, scf, fci
import numpy as np
import os
import sys

destDir="./opencap/tests/pyscf"
sys_dict = {"geometry":'''H 0.0000000000 0.0000000000 0.3705000000
H 0.0000000000 0.0000000000 -0.3705000000''',
            "basis_file":destDir+"/aug-cc-pvdz.bas",
            "molecule": "inline",
            "bohr_coordinates": "false",
            "cart_bf": ""}
cap_dict = {
            "cap_type": "box",
            "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
            "Radial_precision": "14",
            "angular_points": "110"
}

s = pycap.System(sys_dict)
H_bas = gto.basis.load(destDir+'/pyscf.bas', 'H')
mol = gto.M(
            atom = 'H 0.0000000000 0.0000000000 0.3705000000;           \
            H 0.0000000000 0.0000000000 -0.3705000000',
            basis = {'H': H_bas},
            symmetry=True
            )
mol.build()
s = pycap.System(sys_dict)

def test_overlap():
    pyscf_smat = scf.hf.get_ovlp(mol)
    s.check_overlap_mat(pyscf_smat,"pyscf")

def test_pyscf():
    pc = pycap.Projected_CAP(s,cap_dict,5,"pyscf")
    myhf = scf.RHF(mol)
    myhf.kernel()
    fs = fci.FCI(mol, myhf.mo_coeff)
    fs.wfnsym = 'A1g'
    fs.nroots = 5
    fs.spin = 0
    fs.davidson_only = True
    e, c = fs.kernel()
    for i in range(0,len(fs.ci)):
        for j in range(0,len(fs.ci)):
            dm1 = fs.trans_rdm1(fs.ci[i],fs.ci[j],myhf.mo_coeff.shape[1],mol.nelec)
            dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
            pc.add_tdm(dm1_ao,i,j,"pyscf")
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
