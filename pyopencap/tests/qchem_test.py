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
    SOFTWARE.'''

import pyopencap
import numpy as np
import h5py
import os
import sys

destDir="../opencap/tests/qchem"

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

def test_EE():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/EE.fchk"}
    es_dict = {"method" : "EOM",
            "package": "qchem",
        "qchem_output":destDir+"/EE.out",
            "qchem_fchk":destDir+"/EE.fchk",
}
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,cap_dict,3)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_EA():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/EA.fchk"}
    es_dict = {"method" : "EOM",
            "package": "qchem",
        "qchem_output":destDir+"/EA.out",
            "qchem_fchk":destDir+"/EA.fchk",
}
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,cap_dict,3)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_EE_cart():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/EE_cart.fchk"}
    es_dict = {"method" : "EOM",
            "package": "qchem",
        "qchem_output":destDir+"/EE_cart.out",
            "qchem_fchk":destDir+"/EE_cart.fchk",
}
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,cap_dict,3)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_tddft():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/tddft.fchk"}
    es_dict = {"method" : "tddft",
        "package": "qchem",
        "qchem_output":destDir+"/tddft.out",
            "qchem_fchk":destDir+"/tddft.fchk",
}
    s = pyopencap.System(sys_dict)
    pc = pyopencap.CAP(s,cap_dict,4)
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()





