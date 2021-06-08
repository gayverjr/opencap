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
import psi4
import numpy as np
import os
import sys

cap_dict = {
    "cap_type": "box",
        "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
            "Radial_precision": "14",
            "angular_points": "110"
}
molden_dict = { "molecule": "molden", "basis_file": 'h2.molden'}

mol = psi4.geometry("""
    H 0.0000000000 0.0000000000 0.3705000000
    H 0.0000000000 0.0000000000 -0.3705000000""")
E, wfn = psi4.energy('scf/cc-pvdz', return_wfn=True)

mints = psi4.core.MintsHelper(wfn.basisset())
S_mat = np.asarray(mints.ao_overlap())
nstates = 3
psi4.set_options({"opdm":True,"num_roots":nstates,"tdm":True,"dipmom":True})
ci_energy, ci_wfn = psi4.energy('FCI/cc-pvdz', return_wfn=True)
n_bas = S_mat.shape[0]
so2ao = mints.petite_list().sotoao()

def write_molden():
    psi4.molden(wfn, 'h2.molden')
    with open('h2.molden', "a") as myfile:
        myfile.write("\n [7F] \n")

def test_from_molden():
    # add 7F to molden file, psi4 doesn't write it for some reason
    write_molden()
    s = pyopencap.System(molden_dict)
    s.check_overlap_mat(S_mat,"psi4")
    os.remove('h2.molden')

def test_psi4():
    write_molden()
    s = pyopencap.System(molden_dict)
    pc = pyopencap.CAP(s,cap_dict,nstates)
    mo_coeff = ci_wfn.Ca()
    for i in range(0,nstates):
        for j in range(i,nstates):
            opdm_mo = ci_wfn.get_opdm(i, j, "SUM", True)
            opdm_so = psi4.core.triplet(ci_wfn.Ca(), opdm_mo, ci_wfn.Ca(), False, False, True)
            opdm_ao = psi4.core.Matrix(n_bas,n_bas)
            opdm_ao.remove_symmetry(opdm_so,so2ao)
            pc.add_tdm(opdm_ao.to_array(),i,j,"psi4")
            if not i==j:
                pc.add_tdm(opdm_ao.to_array(),j,i,"psi4")
    pc.compute_projected_cap()
    os.remove('h2.molden')
