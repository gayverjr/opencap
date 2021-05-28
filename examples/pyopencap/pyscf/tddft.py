'''Copyright (c) 2020 James Gayvert

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
CAP/TDA-TDDFT(omegaB97XD) on N2-. Basis set is aug-cc-pvtz + 3s3p3d[Gh].
'''

import adcc
import numpy as np
from pyscf import gto,scf,tools,tdscf
import pyopencap
from pandas import DataFrame
import time

start = time.time()

def tda_density_matrix(td, state1, state2):
    '''
        Taking the TDA amplitudes as the CIS coefficients, calculate the density
        matrix (in AO basis) of the excited states
        '''
    ## alpha first
    cis_t1 = td.xy[state1][0][0]
    cis_t2 = td.xy[state2][0][0]
    dm_oo =-np.einsum('ia,ka->ik', cis_t2.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t2.conj())
    
    mf = td._scf
    if state1 == state2:
        # The ground state density matrix in mo_basis
        dm = np.diag(mf.mo_occ[0])
    else:
        dm = np.zeros(mf.mo_coeff[0].shape)
    
    # Add CIS contribution
    nocc = cis_t1.shape[0]
    # Note that dm_oo and dm_vv correspond to spin-up contribution. "*2" to
    # include the spin-down contribution
    dm[:nocc,:nocc] += dm_oo
    dm[nocc:,nocc:] += dm_vv
    
    # Transform density matrix to AO basis
    mo = mf.mo_coeff[0]
    dm_alpha = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    
    ## beta second
    cis_t1 = td.xy[state1][0][1]
    cis_t2 = td.xy[state2][0][1]
    dm_oo =-np.einsum('ia,ka->ik', cis_t2.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t2.conj())
    
    mf = td._scf
    if state1 == state2:
        # The ground state density matrix in mo_basis
        dm = np.diag(mf.mo_occ[1])
    else:
        dm = np.zeros(mf.mo_coeff[1].shape)
    
    # Add CIS contribution
    nocc = cis_t1.shape[0]
    # Note that dm_oo and dm_vv correspond to spin-up contribution. "*2" to
    # include the spin-down contribution
    dm[:nocc,:nocc] += dm_oo * 2
    dm[nocc:,nocc:] += dm_vv * 2
    
    # Transform density matrix to AO basis
    mo = mf.mo_coeff[1]
    dm_beta = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    return dm_alpha,dm_beta

start = time.time()

# define system in pyscf and in pyopencap
ghost_bas = gto.basis.load('n2.nw', 'X')
N_bas = gto.basis.load('n2.nw', 'N')
mol = gto.M(
            atom = 'N        0.0     0.0     0.549; \
            N        0.0     0.0     -0.549;        \
            ghost      0.0     0.0     0.0',
            basis = {'N': N_bas, 'ghost': ghost_bas} )
mol.build()
mf = scf.RHF(mol)
mf.kernel()
pyscf_smat = scf.hf.get_ovlp(mol)
molden_dict = {"basis_file":"molden_in.molden",
    "molecule": "molden"}
tools.molden.from_scf(mf,"molden_in.molden")
s = pyopencap.System(molden_dict)
s.check_overlap_mat(pyscf_smat,"pyscf")

cap_dict = {
    "cap_type": "voronoi",
    "r_cut":"3.00",
}


nstates = 40
pc = pyopencap.CAP(s,cap_dict,nstates)
xc_func = 'WB97XD'

mf = mol.RKS()
mf.xc = xc_func
mf.chkfile = 'n2.chk'
mf.kernel()

# use neutral as initial guess
mol.charge = -1
mol.spin = 1
mol.build()
mf = mol.UKS()
dm = mf.from_chk('n2.chk')
mf.xc = xc_func
mf.run(dm)

mytd = tdscf.TDA(mf).run(nstates=nstates)
print(mytd.e_tot)

H0 = np.zeros((nstates,nstates))
for i in range(0,len(mytd.e_tot)):
    H0[i][i] = mytd.e_tot[i]

for i in range(0,nstates):
    for j in range(i,nstates):
        dm1_alpha,dm1_beta = tda_density_matrix(mytd,i,j)
        pc.add_tdms(dm1_alpha,dm1_beta,i,j,"pyscf")
        if not i==j:
            pc.add_tdms(dm1_alpha,dm1_beta,j,i,"pyscf")

pc.compute_projected_cap()
W=pc.get_projected_cap()

from pyopencap.analysis import CAPHamiltonian as CAPH
my_CAPH = CAPH(H0=H0,W=W)

my_CAPH.export("n2_tddft_opencap.out")

import os
os.remove("molden_in.molden")
os.remove('n2.chk')

end = time.time()
print("Time:")
print(end - start)
