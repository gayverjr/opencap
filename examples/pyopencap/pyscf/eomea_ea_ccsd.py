from pyscf import gto, scf
from pyscf.cc import ccsd
from pyscf.cc.addons import spatial2spin
from pyscf.cc.eom_uccsd import spatial2spin_ea
from pyscf import tools
import numpy as np
import pyopencap
import time
from scipy import linalg as LA

start = time.time()

def partial_trace(rho):
    ''' Trace out the beta spin.

    Parameters
    ------------
    rho: 1RDM in spin-orbital basis

    Returns
    ----------
    rho_traced: spin traced 1RDM

    '''
    nspin = rho.shape[0]
    nbf = int(nspin/2)
    rho_traced = np.zeros((nbf,nbf))
    for i in range(0,nbf):
        for j in range(0,nbf):
            rho_traced[i][j] = np.trace(rho[i*2:i*2+2,j*2:j*2+2])
    return rho_traced

def biorthogonalize(Leigvc,Reigvc,energies):
    ''' Biorthogonalizes a set of left and right eigenvectors from EOM-CC. Each degenerate block is rotated to be biorthogonal separately, and then the entire 
        set is biorthogonalized.

    Parameters
    -----------
    Leigvc: 2d array
        Left eigenvectors as column vectors (i.e. Leigvc.shape[1]=nstates)
    Reigvc: 2d array
        Right eigenvectors as column vectors (i.e. Reigvc.shape[1]=nstates)
    energies: 1d array
        List of eigenenergies
    
    Returns
    --------
    Leigvc: 2d array
        Biorthogonalized left eigenvectors as column vectors (i.e. Leigvc.shape[1]=nstates)
    Reigvc: 2d array
        Biorthogonalized right eigenvectors as column vectors (i.e. Reigvc.shape[1]=nstates) 
    '''
    degen_thresh = 1E-5
    degen_blocks = []
    cur_block = [0]
    cur_energy = energies[0]
    final_Leigvc = np.zeros(Leigvc.shape)
    final_Reigvc = np.zeros(Reigvc.shape)
    for i in range(1,len(energies)):
        if np.abs(energies[i]-cur_energy)<degen_thresh:
            cur_block.append(i)
        else:
            degen_blocks.append(cur_block.copy())
            cur_block = [i]
            cur_energy = energies[i]
    degen_blocks.append(cur_block.copy())
    for block in degen_blocks:
        if len(block) == 1: 
            index = block[0]
            final_Reigvc[:,index] = Reigvc[:,index]
            final_Leigvc[:,index] = Leigvc[:,index]
        else:
            temp_Leigvc = np.zeros((Leigvc.shape[0],len(block)))
            temp_Reigvc = np.zeros((Reigvc.shape[0],len(block)))
            for i,index in enumerate(block):
                temp_Leigvc[:,i] = Leigvc[:,index]
                temp_Reigvc[:,i] = Reigvc[:,index]
            temp_Leigvc, temp_Reigvc = _biorthogonalize(temp_Leigvc,temp_Reigvc)
            for i,index in enumerate(block):
                final_Reigvc[:,index] = temp_Reigvc[:,i]
                final_Leigvc[:,index] = temp_Leigvc[:,i]
    return _biorthogonalize(final_Leigvc,final_Reigvc)


def _biorthogonalize(Leigvc, Reigvc):
    ''' Biorthogonalizes a set of left and right eigenvectors using the LU decomposition. 

    In terms of column vectors:
    M = LU
    Leigvc_new.T = L^inv @ Leigvc.T
    Reigvc_new = Reigvc@ R^inv

    There is an additional check at the end to re-order the left eigenvectors to ensure that the 
    final overlap matrix is unity.

    Parameters
    -----------
    Leigvc: 2d array
        Left eigenvectors as column vectors (i.e. Leigvc.shape[1]=nstates)
    Reigvc: 2d array
        Right eigenvectors as column vectors (i.e. Reigvc.shape[1]=nstates)
    energies: 1d array
        List of eigenenergies
    
    Returns
    --------
    Leigvc: 2d array
        Biorthogonalized left eigenvectors as column vectors (i.e. Leigvc.shape[1]=nstates)
    Reigvc: 2d array
        Biorthogonalized right eigenvectors as column vectors (i.e. Reigvc.shape[1]=nstates) 
    '''
    M = Leigvc.T @ Reigvc
    P, L, U = LA.lu(M)
    Linv = LA.inv(L)
    Uinv = LA.inv(U)
    Leigvc = np.dot(Linv, Leigvc.T)
    Reigvc = np.dot(Reigvc, Uinv)
    Leigvc = Leigvc.T
    # re-order left eigenvectors just in case
    M = Leigvc.T @ Reigvc
    L_reordered = np.zeros(Leigvc.shape)
    for i in range(0,len(M)):
        jmax = np.argmax(Leigvc.T@Reigvc[:,i])
        L_reordered[:,i] = Leigvc[:,jmax]
    return L_reordered, Reigvc

def make_rdm1(t1,t2,r1,r2,l1,l2,state_dm=False):
    d1 = _gamma1_intermediates(t2,r1,r2,l2)
    return _make_rdm1(mycc,d1,t1,t2,r1,r2,l1,state_dm)


def _gamma1_intermediates(t2,r1,r2,l2):
    lt_oo = 0.5 * np.einsum('jcd,icd->ji',l2,r2)
    lt_vv = np.einsum('kbc,kac->ba',l2,r2)
    y1_ov = np.einsum('iac,c->ia',l2,r1)
    ltt_o = 0.5 * np.einsum('kdc,ikdc->i',l2,t2)
    return lt_oo,lt_vv,y1_ov,ltt_o


def _make_rdm1(d1,t1,t2,r1,r2,l1,state_dm=False):
    lt_oo,lt_vv,y1_ov,ltt_o = d1
    nocc,nvir = y1_ov.shape
    nmo = nocc + nvir
    dm1 = np.zeros((nmo,nmo),dtype=lt_oo.dtype)
    if state_dm:
        dm1[np.diag_indices(nocc)] += 1
    ## OO
    dm1[:nocc,:nocc] -= np.einsum('ji->ij',lt_oo)
    dm1[:nocc,:nocc] -= np.einsum('jc,ic->ij',y1_ov,t1)
    ## OV
    if state_dm:
        dm1[:nocc,nocc:] += np.einsum('ia->ia',t1)
    dm1[:nocc,nocc:] -= np.einsum('ic,ca->ia',t1,lt_vv)
    dm1[:nocc,nocc:] -= np.einsum('ka,ki->ia',t1,lt_oo)
    dm1[:nocc,nocc:] -= np.einsum('i,a->ia',ltt_o,r1)
    dm1[:nocc,nocc:] += np.einsum('iac,a->ia',r2-np.einsum('ic,a->ica',t1,r1),l1)
    dm1[:nocc,nocc:] += np.einsum('kica,kc->ia',t2-np.einsum('ka,ic->kica',t1,t1),y1_ov)
    ## VV
    dm1[nocc:,nocc:] += np.einsum('b,a->ab',r1,l1)
    dm1[nocc:,nocc:] += np.einsum('ab->ab',lt_vv)
    dm1[nocc:,nocc:] += np.einsum('ka,kb->ab',y1_ov,t1)
    ## VO
    dm1[nocc:,:nocc] -= np.einsum('ia->ai',y1_ov)
    return dm1 
 

ghost_bas = gto.basis.load('n2.nw', 'X')
N_bas = gto.basis.load('n2.nw', 'N')
mol = gto.M(
            atom = '    N 0.0000000000 0.0000000000 0.548756750;\
            N 0.0000000000 0.0000000000 -0.548756750;\
            ghost      0.0     0.0     0.0',
            basis = {'N': N_bas, 'ghost': ghost_bas} )
mol.verbose = 1
mol.build()
mf = mol.RHF().run()
# create system using molden file
molden_dict = {"basis_file":"molden_in.molden",
    "molecule": "molden"}
tools.molden.from_scf(mf,"molden_in.molden")
s = pyopencap.System(molden_dict)
pyscf_smat = scf.hf.get_ovlp(mol)
s.check_overlap_mat(pyscf_smat,"pyscf")
cap_dict = {"cap_type": "box",
    "cap_x":"2.76",
    "cap_y":"2.76",
    "cap_z":"4.88",
}


nstates = 32

pc = pyopencap.CAP(s,cap_dict,nstates)
pc.compute_ao_cap(cap_dict)


mycc = ccsd.CCSD(mf)
mycc.kernel()
t1 = mycc.t1
t2 = mycc.t2
nocc, nvir = t1.shape
# convert to spin orbital basis
t1_spin = spatial2spin(t1)
t2_spin = spatial2spin(t2)

myeom = mycc.eomea_method()
e_r, r_eom = mycc.eaccsd(nroots=nstates)
e_l, l_eom = mycc.eaccsd(nroots=nstates,left=True)
# use as column vectors
r_eom = np.array(r_eom).T
l_eom = np.array(l_eom).T
l_eom,r_eom = biorthogonalize(l_eom,r_eom,e_r)

r1 = []
r2 = []
l1 = []
l2 = []
h0 = np.zeros((nstates,nstates))
# collect EOM amplitudes in spin orbital basis
for i in range(0,nstates):
    r1_curr,r2_curr = myeom.vector_to_amplitudes(r_eom[:,i],myeom.nmo,myeom.nocc)
    l1_curr,l2_curr = myeom.vector_to_amplitudes(l_eom[:,i],myeom.nmo,myeom.nocc)
    dummyR1 = np.zeros(np.shape(r1_curr))
    dummyR2 = np.zeros(np.shape(r2_curr))
    r1_spin,r2_spin = spatial2spin_ea((r1_curr,dummyR1),(dummyR2,r2_curr,dummyR2,dummyR2))
    l1_spin,l2_spin = spatial2spin_ea((l1_curr,dummyR1),(dummyR2,l2_curr,dummyR2,dummyR2))
    r1.append(r1_spin)
    r2.append(r2_spin)
    l1.append(l1_spin)
    l2.append(l2_spin)
    h0[i][i] = e_r[i] + mycc.e_hf + mycc.e_corr

for i in range(0,nstates):
    for j in range(0,nstates):
        if i == j:
            dm1 = make_rdm1(mycc,t1_spin,t2_spin,r1[i],r2[i],l1[i],l2[i],state_dm=True)	
        else:
            dm1 = make_rdm1(mycc,t1_spin,t2_spin,r1[i],r2[i],l1[j],l2[j],state_dm=False)
        dm1 = partial_trace(dm1)
        dm1_ao = mf.mo_coeff @ dm1 @ mf.mo_coeff.T
        pc.add_tdm(dm1_ao,i,j,'pyscf')

pc.compute_projected_cap()
W=pc.get_projected_cap()

from pyopencap.analysis import CAPHamiltonian as CAPH
my_CAPH = CAPH(H0=h0,W=W)

my_CAPH.export("n2_eomea.out")

end = time.time()
print("Time:")
print(end - start)
