import pyopencap as pycap
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
    es_dict = {"method" : "EOMEE",
        "qchem_output":destDir+"/EE.out",
            "qchem_fchk":destDir+"/EE.fchk",
}
    s = pycap.System(sys_dict)
    pc = pycap.Projected_CAP(s,cap_dict,3,"qchem")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_EA():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/EA.fchk"}
    es_dict = {"method" : "EOMEA",
        "qchem_output":destDir+"/EA.out",
            "qchem_fchk":destDir+"/EA.fchk",
}
    s = pycap.System(sys_dict)
    pc = pycap.Projected_CAP(s,cap_dict,3,"qchem")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_IP():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/IP.fchk"}
    es_dict = {"method" : "EOMIP",
        "qchem_output":destDir+"/IP.out",
            "qchem_fchk":destDir+"/IP.fchk",
}
    s = pycap.System(sys_dict)
    pc = pycap.Projected_CAP(s,cap_dict,3,"qchem")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()

def test_EE_cart():
    sys_dict = {"molecule": "qchem_fchk",
"basis_file": destDir+"/EE_cart.fchk"}
    es_dict = {"method" : "EOMEE",
        "qchem_output":destDir+"/EE_cart.out",
            "qchem_fchk":destDir+"/EE_cart.fchk",
}
    s = pycap.System(sys_dict)
    pc = pycap.Projected_CAP(s,cap_dict,3,"qchem")
    pc.read_data(es_dict)
    pc.compute_ao_cap()
    pc.renormalize()
    pc.compute_projected_cap()
    mat=pc.get_projected_cap()
    h0 = pc.get_H()






