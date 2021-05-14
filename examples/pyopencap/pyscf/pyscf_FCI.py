import pyopencap
from pandas import DataFrame
from pyscf import gto, scf, ci, ao2mo, tools, fci
import numpy as np

cap_dict = {
            "cap_type": "box",
            "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
}
nstates = 22

# pyscf mol
ghost_bas = gto.basis.load('H2_custom.bas', 'X')
H_bas = gto.basis.load('H2_custom.bas', 'H')
mol = gto.M(
            atom = 'ghost 0.0000000000 0.0000000000 -0.7400000000; \
            ghost 0.0000000000 0.0000000000 0.7400000000;          \
            H 0.0000000000 0.0000000000 0.3705000000;           \
             H 0.0000000000 0.0000000000 -0.3705000000',
             basis = {'H': H_bas, 'ghost': ghost_bas},
             symmetry=True)
mol.build()
myhf = scf.RHF(mol)
myhf.kernel()

molden_dict = {"basis_file":"molden_in.molden",
    "molecule": "molden"}
tools.molden.from_scf(myhf,"molden_in.molden")
s = pyopencap.System(molden_dict)
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")

# run full ci, takes a while
fs = fci.FCI(mol, myhf.mo_coeff)
fs.nroots = nstates
fs.spin = 0
e, c = fs.kernel()
print("Ground state energy:" + str(e[0]))

# create cap object
pc = pyopencap.CAP(s,cap_dict,nstates,"pyscf")
pc.compute_ao_cap()

# fill density matrices
h0 = np.zeros((nstates,nstates))
for i in range(0,len(e)):
    h0[i][i]=e[i]
for i in range(0,len(fs.ci)):
    for j in range(0,len(fs.ci)):
        dm1 = fs.trans_rdm1(fs.ci[i],fs.ci[j],myhf.mo_coeff.shape[1],mol.nelec)
        dm1_ao = np.einsum('pi,ij,qj->pq', myhf.mo_coeff, dm1, myhf.mo_coeff.conj())
        pc.add_tdm(dm1_ao,i,j,"pyscf")

pc.compute_projected_cap()
mat=pc.get_projected_cap()

print("Printing out matrices required for Projected CAP calculation.")
print("Number of states: "+ str(nstates))
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP matrix")
print(DataFrame(mat).to_string(index=False, header=False))

