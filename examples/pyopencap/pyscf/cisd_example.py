#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
"""Example computing the TPA cross section for water using ADC
    (10.1063/1.3682324)
    """
import numpy as np
from pyscf import gto,scf,tools,dft,tdscf,mcscf,tddft,ci,fci
import pyopencap
import pyscf
from pandas import DataFrame

cap_dict = {
    "cap_type": "voronoi",
        "r_cut": "5.00",
            "Radial_precision": "18",
            "angular_points": "590"
}


mol=tools.molden.load('be.molden')[0]
mol.charge = 0
mol.spin = 0
mol.verbose=4
mol.build()

mf = scf.RHF(mol)
mf.kernel()
mycc = mf.CISD()
mycc.kernel()


nstates = 15
mol=tools.molden.load('be.molden')[0]
mol.charge = -1
mol.spin = 1
mol.verbose=4
mol.build()

mf = scf.UHF(mol)
mf.kernel()
mycc = mf.CISD()
mycc.nstates = nstates
mycc.kernel()




nstates = len(mycc.e_tot)

molden_dict = {"basis_file":"be.molden",
    "molecule": "molden"}
s = pyopencap.System(molden_dict)
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")
pc = pyopencap.CAP(s,cap_dict,nstates,"pyscf")

h0 = np.zeros((nstates,nstates))
for i in range(0,nstates):
    h0[i][i]= mycc.e_tot[i]

print(mycc.e_tot)

for i in range(0,nstates):
    for j in range(0,nstates):
        dm_a, dm_b = mycc.trans_rdm1(mycc.ci[i],mycc.ci[j])
        dma_ao = np.einsum('pi,ij,qj->pq', mf.mo_coeff[0], dm_a, mf.mo_coeff[0].conj())
        dmb_ao = np.einsum('pi,ij,qj->pq', mf.mo_coeff[1], dm_b, mf.mo_coeff[1].conj())
        pc.add_tdms(dma_ao,dmb_ao,i,j,"pyscf")


pc.compute_ao_cap()
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()


print("Printing out matrices required for Perturbative CAP calculation.")
print("Number of states: " + str(nstates))
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP Matrix")
print(DataFrame(mat).to_string(index=False, header=False))




