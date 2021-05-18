import psi4
from pandas import DataFrame
import numpy as np
import pyopencap
import time

start = time.time()
psi4.core.set_output_file('output.dat', False)

mol = psi4.geometry("""
    @He 0.0000000000 0.0000000000 -0.7400000000
    @He 0.0000000000 0.0000000000 0.7400000000
    H 0.0000000000 0.0000000000 0.3705000000
    H 0.0000000000 0.0000000000 -0.3705000000
    Symmetry C1""")
def basisspec_psi4_yo__anonymous03952cbd(mol, role):
    basstrings = {}
    mol.set_basis_all_atoms("test_bas", role=role)
    basstrings['test_bas'] = """
        spherical
        ****
        H    0
        S   3  1.00
        82.6400000000          0.00200600
        12.4100000000          0.01534300
        2.8240000000           0.07557900
        S   1  1.00
        0.7977000000           1.00000000
        S   1  1.00
        0.2581000000           1.00000000
        S   1  1.00
        0.0898900000           1.00000000
        S   1  1.00
        0.0236300000           1.00000000
        P   1  1.00
        2.2920000000           1.00000000
        P   1  1.00
        0.8380000000           1.00000000
        P   1  1.00
        0.2920000000           1.00000000
        P   1  1.00
        0.0848000000           1.00000000
        D   1  1.00
        2.0620000000           1.00000000
        D   1  1.00
        0.6620000000           1.00000000
        D   1  1.00
        0.1900000000           1.00000000
        F   1  1.00
        1.3970000000           1.00000000
        F   1  1.00
        0.3600000000           1.00000000
        ****
        He 0
        S 1  1.00
        0.0118150000 1.00000000
        S 1  1.00
        0.0059075000 1.00000000
        S 1  1.00
        0.0029537500 1.00000000
        S 1  1.00
        0.0014768750 1.00000000
        P 1  1.00
        0.0424000000 1.00000000
        P 1  1.00
        0.0212000000 1.00000000
        P 1  1.00
        0.0106000000 1.00000000
        P 1  1.00
        0.0053000000 1.00000000
        ****
        """
    return basstrings

nstates = 10
psi4.qcdb.libmintsbasisset.basishorde['ANONYMOUS03952CBD'] = basisspec_psi4_yo__anonymous03952cbd
psi4.core.set_global_option("BASIS", "anonymous03952cbd")
E, wfn = psi4.energy('scf', return_wfn=True)
mints = psi4.core.MintsHelper(wfn.basisset())
S_mat = np.asarray(mints.ao_overlap())
n_bas = S_mat.shape[0]
so2ao = mints.petite_list().sotoao()
psi4.molden(wfn, 'h2.molden')
# add 7F to molden file, psi4 doesn't write it for some reason
with open("h2.molden", "a") as myfile:
    myfile.write("\n [7F] \n")

molden_dict = {"basis_file":"h2.molden",
    "molecule": "molden"}
s = pyopencap.System(molden_dict)
s.check_overlap_mat(S_mat,"psi4")
cap_dict = {
    "cap_type": "box",
        "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
}
pc = pyopencap.CAP(s,cap_dict,nstates)
pc.compute_ao_cap()

psi4.set_options({"opdm":True,"num_roots":nstates,"tdm":True,"dipmom":True})
ci_energy, ci_wfn = psi4.energy('FCI', return_wfn=True)
h0 = np.zeros((nstates,nstates))
for i in range(0,nstates):
    var_str = 'CI ROOT ' + str(i) + ' TOTAL ENERGY'
    h0[i][i] = float(psi4.get_variable(var_str))

mo_coeff = ci_wfn.Ca()
for i in range(0,nstates):
    for j in range(i,nstates):
        opdm_mo = ci_wfn.get_opdm(i, j, "SUM", True)
        opdm_ao = psi4.core.triplet(mo_coeff, opdm_mo, mo_coeff, False, False, True)
        pc.add_tdm(opdm_ao.to_array(),i,j,"psi4")
        if not i==j:
            pc.add_tdm(opdm_ao.to_array(),j,i,"psi4")

pc.compute_projected_cap()
mat=pc.get_projected_cap()

print("Printing out matrices required for Projected CAP calculation.")
print("Number of states: "+ str(nstates))
print("Zeroth order Hamiltonian")
print(DataFrame(h0).to_string(index=False, header=False))
print("CAP matrix")
print(DataFrame(mat).to_string(index=False, header=False))

end = time.time()
print("Time:")
print(end-start)

