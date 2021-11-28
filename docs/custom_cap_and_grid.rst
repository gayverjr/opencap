.. _custom:

Custom CAPs and Grids
=======================================
Starting with PyOpenCAP version 1.2, users can now specify customized CAP functions and numerical integration grids. 

Custom CAPs
------------
Python functions with the signature ``vector<double>,vector<double>,vector<double>,vector<double> --> vector<double>`` 
can be used as CAP functions by the :class:`~pyopencap.CAP` class for numerical integration. An example is provided below:

.. code-block:: python

    # this defines a box CAP of with cutoffs of 3 bohr in each coordinate
    def box_cap(x,y,z,w):
        cap_values = []
        cap_x = 3.00
        cap_y = 3.00
        cap_z = 3.00
        for i in range(0,len(x)):
            result = 0
            if np.abs(x[i])>cap_x:
                result += (np.abs(x[i])-cap_x) * (np.abs(x[i])-cap_x)
            if np.abs(y[i])>cap_y:
                result += (np.abs(y[i])-cap_y) * (np.abs(y[i])-cap_y)
            if np.abs(z[i])>cap_z:
                result += (np.abs(z[i])-cap_z) * (np.abs(z[i])-cap_z)
            result = w[i]*result
            cap_values.append(result)
        return cap_values

    cap_dict = {"cap_type": "custom"}
    pc = pyopencap.CAP(s,cap_dict,5,box_cap)

Custom Grids
--------------
Custom grids for numerical integration can be specified using the :class:`~pyopencap.CAP.compute_cap_on_grid` function. The arguments are 
assumed to be 1D arrays of equal size. The function can be called repeatedly for a cumulative sum in the case of atomic grids. An example is provided below:

.. code-block:: python

    for i in range(0,Natoms):
        x,y,z,w = get_grid_for_atom(atoms[i])
        pc.compute_cap_on_grid(x,y,z,w)
    # final sum is cumulative
    pc.compute_projected_cap()

