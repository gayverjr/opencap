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
import numpy as np
import os
import sys
import pytest

try:
    import psi4
except: 
    pass


cap_dict = {
    "cap_type": "box",
        "cap_x":"6.00",
            "cap_y":"6.00",
            "cap_z":"6.7",
            "Radial_precision": "14",
            "angular_points": "110"
}
molden_dict = { "molecule": "molden", "basis_file": 'h2.molden'}
nstates = 3

def get_ci_wfn():
    mol = psi4.geometry("""
        H 0.0000000000 0.0000000000 0.3705000000
        H 0.0000000000 0.0000000000 -0.3705000000""")
    E, wfn = psi4.energy('scf/cc-pvdz', return_wfn=True)
    psi4.set_options({"opdm":True,"num_roots":nstates,"tdm":True,"dipmom":True})
    ci_energy, ci_wfn = psi4.energy('FCI/cc-pvdz', return_wfn=True)
    mints = psi4.core.MintsHelper(wfn.basisset())
    return ci_wfn, mints
    
def write_molden(wfn):
    psi4.molden(wfn, 'h2.molden')
    # add 7F to molden file, psi4 doesn't write it for some reason
    with open('h2.molden', "a") as myfile:
        myfile.write("\n [7F] \n")

@pytest.mark.skipif('psi4' not in sys.modules,reason="requires the Psi4 library")
def test_from_molden():
    wfn, mints = get_ci_wfn()
    write_molden(wfn)
    S_mat = np.asarray(mints.ao_overlap())
    s = pyopencap.System(molden_dict)
    s.check_overlap_mat(S_mat,"psi4")
    os.remove('h2.molden')

@pytest.mark.skipif('psi4' not in sys.modules,reason="requires the Psi4 library")
def test_psi4():
    ci_wfn, mints = get_ci_wfn()
    write_molden(ci_wfn)
    s = pyopencap.System(molden_dict)
    pc = pyopencap.CAP(s,cap_dict,nstates)
    mo_coeff = ci_wfn.Ca()
    S_mat = np.asarray(mints.ao_overlap())
    n_bas = S_mat.shape[0]
    so2ao = mints.petite_list().sotoao()
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
