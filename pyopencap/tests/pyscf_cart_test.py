import pyopencap
from pyscf import gto, scf, fci
import numpy as np
import os
import sys

destDir="../opencap/tests/data"
sys_dict = {"geometry":'''H        0.0     0.0     0.54981512
    H        0.0     0.0     -0.54981512
    X       0.0     0.0     0.0''',
        "basis_file":destDir+"/test_basis_for_pyscf.bas",
        "molecule": "inline",
        "cart_bf" : "dfg"
}
cap_dict = {
    "cap_type": "box",
        "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
            "Radial_precision": "14",
            "angular_points": "110"
}

s = pyopencap.System(sys_dict)
H_bas = gto.basis.load(destDir+'/pyscf_test_basis.bas', 'H')
X_bas = gto.basis.load(destDir+'/pyscf_test_basis.bas', 'O')
mol = gto.M(
            atom = 'H        0.0     0.0     0.54981512;\
            H        0.0     0.0     -0.54981512;\
            X       0.0     0.0     0.0',
            basis = {'H': H_bas, 'X':X_bas}, cart=True
            )
mol.build()

def test_overlap():
    basis_id_str = s.get_basis_ids()
    pyscf_smat = scf.hf.get_ovlp(mol)
    s.check_overlap_mat(pyscf_smat,"pyscf")

def test_pyscf():
    pc = pyopencap.Projected_CAP(s,cap_dict,3,"pyscf")
    myhf = scf.RHF(mol)
    myhf.kernel()
    fs = fci.FCI(mol, myhf.mo_coeff)
    fs.nroots = 3
    e, c = fs.kernel()
    for i in range(0,len(fs.ci)):
        for j in range(0,len(fs.ci)):
            dm1 = fs.trans_rdm1(fs.ci[i],fs.ci[j],myhf.mo_coeff.shape[1],mol.nelec)
            dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
            pc.add_tdm(dm1_ao,i,j,"pyscf")
    pc.compute_ao_cap()
    pyscf_smat = scf.hf.get_ovlp(mol)
    pc.renormalize_cap(pyscf_smat,"pyscf")
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
