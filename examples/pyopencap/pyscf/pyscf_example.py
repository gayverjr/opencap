import pyopencap
from pandas import DataFrame
from pyscf import gto, scf, ci, ao2mo, tools, fci
import numpy as np
sys_dict = {"geometry":    '''X 0.0000000000 0.0000000000 -0.7400000000
X 0.0000000000 0.0000000000 0.7400000000
H 0.0000000000 0.0000000000 0.3705000000
H 0.0000000000 0.0000000000 -0.3705000000''',
            "basis_file":"psi4.bas",
            "molecule": "inline"}
cap_dict = {
            "cap_type": "box",
            "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
            "Radial_precision": "14",
            "angular_points": "110"
}
nstates = 22

# pyscf mol
ghost_bas = gto.basis.load('nwchem.bas', 'N')
H_bas = gto.basis.load('nwchem.bas', 'H')
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

# create system, check overlap matrix
s = pyopencap.System(sys_dict)
my_smat = s.get_overlap_mat()
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")


# we can also create the system using a molden file
molden_dict = {"basis_file":"molden_in.molden",
    "molecule": "molden"}
tools.molden.from_scf(myhf,"molden_in.molden")
s = pyopencap.System(molden_dict)
s.check_overlap_mat(pyscf_smat,"pyscf")


# run full ci, takes a while
fs = fci.FCI(mol, myhf.mo_coeff)
fs.wfnsym = 'A1g'
fs.nroots = nstates
fs.spin = 0
e, c = fs.kernel()
print("Ground state energy:" + str(e[0]))

# create projected cap object
pc = pyopencap.Projected_CAP(s,cap_dict,nstates,"pyscf")
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

# select the states we want
save_indices = [2,5,6,8,10,12,14,16,18,20,22]
for i in range(22,-1,-1):
    if i not in save_indices:
        h0 = np.delete(h0,i,axis=0)
        h0 = np.delete(h0,i,axis=1)
        mat = np.delete(mat,i,axis=0)
        mat = np.delete(mat,i,axis=1)
print("Printing out matrices required for Projected CAP calculation.")
print("Number of states: 10")
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP matrix")
print(DataFrame(mat).to_string(index=False, header=False))
