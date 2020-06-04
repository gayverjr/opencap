import numpy as np
from pyscf import gto, scf, ci, ao2mo
import pycap
from pandas import DataFrame


sys_dict = {"geometry":    '''N  0  0   1.039
    N  0  0   -1.039
    Gh 0  0   0.0''',
        "basis_file":"test_bas.bas",
            "bohr_coordinates": "true",
            "cart_bf": ""}
sys_dict2 = {"geometry":    '''N  0  0   0.55
    N  0  0   -0.55''',
            "basis_file":"test_bas.bas",
            "bohr_coordinates": "false",
            "cart_bf": ""}

# pyscf stuff
s = pycap.System(sys_dict2)
Smat1 = s.get_overlap_mat()
#print(DataFrame(Smat1).to_string(index=False, header=False))

ghost_bas = gto.basis.load('/Users/JG/Documents/Software/OpenCAP/opencap/example/pyscf_bas.bas', 'H')
g2= gto.basis.load('/Users/JG/Documents/Software/OpenCAP/opencap/example/pyscf_bas.bas', 'N')



mol = gto.M(
            atom = 'N 0 0 0.55; N 0 0 -0.55; ghost 0 0 0',
            basis = {'N': g2, 'ghost': ghost_bas}
            )

mol = gto.M(
            atom = 'N 0 0 0.55; N 0 0 -0.55',
            basis = {'N': g2, 'ghost': ghost_bas}
            )

Smat2 = scf.hf.get_ovlp(mol)
#print(DataFrame(Smat2).to_string(index=False, header=False))

for i in range(0,len(Smat1)):
    for j in range(0,len(Smat1)):
        if np.abs(Smat1[i][j]-Smat2[i][j])>=1E-2:
            print(str(i)+"," + str(j)+":" + str(Smat1[i][j]-Smat2[i][j]))
            print("Smat1:" + str(Smat1[i][j]))
            print("Smat2:" + str(Smat2[i][j]))
            print("\n")



'''
mf = scf.RHF(mol).run()
myci = ci.CISD(mf).run()
myci.nroots = 6
myci.kernel()

t_dm1 = myci.trans_rdm1(myci.ci[3], myci.ci[0])
t_dm1 = numpy.einsum('pi,ij,qj->pq', mf.mo_coeff, t_dm1, mf.mo_coeff.conj())
'''

