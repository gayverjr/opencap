Input {#input}
==============

Input sections
==============

An OpenCAP input file has four sections: 

$system and $projected_cap are required, $geometry and $trajectory are optional.
  
Within each section, the user inserts keyword/value pairs which define the calculation.  
All four sections have their own unique keywords, and all four must be specified to run a calculation. 
The general format is:

    $section_name

    keyword  value

    keyword  value

    $end

In each section, the parser will read lines until it sees the $end keyword.  
Only the first two arguments per line are read, so all remaining text can be used as comments. 
Additionally, lines beginning with an exclamation mark '!' are ignored. All input sections 
(except for those referring to file paths) are case-insensitive.

System
=======

The system section specifies information about the molecular system and the ab initio basis set. Basis sets can be downloaded from the [MolSSI basis set exchange](https://www.basissetexchange.org/) and modified to suit your purposes.

__Required__
| Keyword    | Valid options                     | Description                                                                                                                                                                                                                                                                                      |
|------------|-----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| molecule   | molden,qchem_fchk,rassi_h5,inline | Specifies which format to read the molecular geometry. If "inline" is chosen, the "$geometry" section is also required.                                                                                                                                                                          |
| basis_file | path to basis file                | Specifies the path to the basis file. When "molecule" is set to "molden","rassi_h5", or "qchem_fchk", this keyword should be set to a path to a file of the specified type. When "molecule" is set to  "inline", this keyword should be set to a path to a basis set file formatted in "Psi4" style. |

__Optional__
| Keyword          | default                | Description                                                                                                                                                                                                                                                                                                                                                                                |
|------------------|------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| cart_bf          | ""(all functions pure) | Controls the use of pure or Cartesian angular forms of GTOs. The letters corresponding to the angular momenta listed in this field will be expanded in cartesians, those not listed will be expanded in pure GTOs. For example, "df" means d and f-type functions will be cartesian, and all others will be pure harmonic. This keyword is only active when "molecule" is set to "inline". |
| bohr_coordinates | false                  | Set to true when coordinates specified in "geometry" keyword are in bohr units.  This keyword is only active when "molecule" is set to "inline".                                                                                                                                                                                                                                           |

__Geometry__

When specifying the geometry inline, use the following format:

    atom1 x-coordinate y-coordinate z-coordinate 

    atom2 x-coordinate y-coordinate z-coordinate ...

Ghost centers with zero nuclear charge can be specified using the symbol "X".

By default, units are assumed to be Angstroms.

Projected_CAP
============

This section allows one to specify the parameters used to define the complex absorbing potential, 
the data read in from the electronic structure calculation, and the grid used for numerical integration.

There are two key pieces of data which must be read in from the electronic structure package: 
the set of one particle transition densities between each pair of states 
(each has a dimension of NxN where N is the number of basis functions), 
and the zeroth order Hamiltonian (which has a dimension of MxM, 
where M is the number of electronic states computed by the package). 
In some cases, the zeroth order Hamiltonian can be read directly from the output file of 
the electronic structure package, but in other cases it may be easier to just specify the 
zeroth order hamiltonian manually (due to its small dimensionality). For the latter 
option, the 'h0_file' keyword must be used, and the value must point to a properly formatted file. 

A properly formatted zeroth order hamiltonian file can have one of two forms. The first is the diagonal form (where the zeroth order hamiltonian just consists of total energies on the diagonal, as is the case for EOM-CCSD), the word "Diagonal" must appear on the first line, followed by the energies, each on their own line. The second is the full form, in which the entire zeroth order hamiltonian is specified. The word "Full" must appear on the first line, followed by the matrix. The matrix is read in row-major order.
We currently support two types of complex absorbing potentials; quadratic box-type potentials, and a quadratic "Voronoi isosurface" potential. Both types of %CAP have their own set of parameters which must be defined. Numerical integration is performed with the aid of the [numgrid](https://github.com/dftlibs/numgrid) package, which allocates the grid based on the molecular geometry and basis set. By default, a radial grid with precision 1e-14 is used, with an angular grid of 590 points. Both can optionally be adjusted to increase or decrease the size of the grid as needed.

__Required__
| Keyword | Type | Description | Valid options |
|:-------------:|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| cap_type | string | Type of absorbing potential| 'Box' or 'Voronoi'|
| Method | string | Name of electronic structure method | MS-CASPT2, XMS-CASPT2, TDDFT, EOM, SC-NEVPT2, PC-NEVPT2|
| Package | string | Name of electronic structure package | OpenMolcas, Q-Chem |
| Nstates | int | Number of electronic states | 2 or greater |

__OpenMolcas Interface__
| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
| rassi_h5 | string | Relative or absolute path to rassi.h5 file containing TDMs.  |
| molcas_output | string | Relative or absolute path to OpenMolcas output file.  Required when h0_file keyword is unspecified.  |

__QChem Interface__
| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
| qchem_fchk | string | Relative or absolute path to qchem .fchk file containing TDMs.  |
| qchem_output | string | Relative or absolute path to Q-Chem output file.  Required when h0_file keyword is unspecified.  |

__Box %CAP__

A quadratic potential which encloses the system in a 3D rectangular box. Functional form:

\f$    W= W_x + W_y +W_z  \f$ 

\f$ W_{\alpha} = 0\f$, \f$ \left|r_{\alpha}\right| < R_{\alpha}^0 \f$

\f$ W_{\alpha} = \left(r_{\alpha} - R_{\alpha}^0 \right)^2,\f$  \f$ \left|r_{\alpha}\right| > R_{\alpha}^0 \f$

| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
| CAP_X | float | Onset of %CAP in x-direction. Specify in bohr units.|
| CAP_Y | float | Onset of %CAP in y-direction. Specify in bohr units.|
| CAP_Z | float | Onset of %CAP in z-direction. Specify in bohr units.|

__Voronoi %CAP__

A quadratic potential which uniformly wraps around the system at a specified cutoff radius. 
The edges between between Voronoi cells are smoothed out to make the potential more 
amenable to numerical integration[1]. Functional form:

\f$ W(\vec{r}) = 0,\f$  \f$r_{WA} \leq r_{cut} \f$ 

\f$ W(\vec{r}) = (r_{WA} - r_{cut} )^2,\f$   \f$r_{WA} > r_{cut} \f$

\f$ r_{WA}(\vec{r}) = \sqrt{\frac{\sum_{i} w_{i}|\vec{r}-\vec{R}_i|^2}{\sum_{i} w_{i}}} \f$
    
\f$ w_{i} = \frac{1}{(|\vec{r}-\vec{R}_i|^2-r_{min}^2+1 a.u.)^2} \f$
	
\f$ r_{min} = \min\limits_{i}{|\vec{r}-\vec{R}_i|} \f$

| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
|r_cut | float | Cutoff radius for Voronoi %CAP. Specify in bohr units.|


__Optional__
| Keyword | Type | Description | Default |
|------------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| radial_precision| int | Radial precision for numerical integration grid. A precision of 1x10^(-N), where N is the value specified is used.  | 14  |
| angular_points | int | Number of angular ponts used for the grid. See  https://github.com/dftlibs/numgrid for allowed numbers of points. | 590 |
| h0_file | string | Relative or absolute path to formatted zeroth order Hamiltonian file.  Required when output from electronic structure package (e.g. molcas_output) is unspecified. | "" |

Trajectory
============
One can generate eigenvalue trajectories within OpenCAP using the $trajectory field. One can 
specify the range of \f$\eta\f$ values to be explored, and for each state, the minimum of 
the logarithmic velocity will be printed in the output. States are tracked using eigenvector overlap.
The entire trajectory can also be exported to a .data file using the save_trajectory keyword.

| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
| eta_step | float | \f$\eta\f$-step size, in 1E-5 bohr|
| nsteps | int | Number of steps for eigenvalue trajectory|
| save_trajectory | bool | Set to true if you want to save the trajectory, it will be saved to a .data file.|

References
==========
[1] Sommerfeld, T.; Ehara, M. Complex Absorbing Potentials with Voronoi Isosurfaces Wrapping Perfectly around Molecules. *J. Chem. Theory Comput.* **2015**, 11 (10), 4627–4633.
