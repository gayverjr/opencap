import pycap
my_dict = {"geometry":    '''N  0  0   1.039
                             N  0  0   -1.039
                            Gh 0  0   0.0''',
            "basis_file":"test_bas.bas",
            "bohr_coordinates": "true",
            "cart_bf": "",
            "method" : "ms-caspt2",
            "package": "openmolcas",
            "molcas_output":"anion_reference.out",
            "rassi_h5":"test.rassi.h5",
            "nstates": "10",
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "radial_precision": "14",
            "angular_points": "110"}
#s = pycap.System()
s = pycap.System(my_dict)
pc = pycap.Projected_CAP(s)
pc.run()