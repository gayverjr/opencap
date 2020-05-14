# OpenCAP input structure #

An OpenCAP input file has four sections; $geometry, $system, $projected_cap, and $job.  Within each section, the user inserts keyword/value pairs which define the calculation.  All four sections have their own unique keywords, and all four must be specified to run a calculation. The general format is:

$section_name

keyword  value

keyword  value

$end

In each section, the parser will read lines until it sees the $end keyword.  Only the first two arguments per line are read, so all remaining text can be used as comments. Additionally, lines beginning with an exclamation mark '!' are ignored. All input sections are case-insensitive.

##  Job ##
This section specifies what type of job is to be executed by the software. We currently only support Projected %CAP calculations, but there's more to come!

#### Required ####
| Keyword | Type | Description | Valid options |
|:-------------:|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| title | string | Title of job | N/A |
| jobtype | string | Type of job to be executed | Projected_cap |

## Geometry ##
Molecular geometries are specified inline as follows:

atom x-coordinate y-coordinate z-coordinate

By default, coordinates are assumed to be in angstroms. To specify in bohr units, one must
set bohr_coordinates to true in the $system section. Ghost atoms (which have 0 nuclear charge but can be used to augment the basis set with additional basis functions functions), are specified in the geometry section with the label "Gh". 

## System ##
The system section specifies information about the molecular system and the ab initio basis set.

#### Required ####
| Keyword | Type | Description |
|:-------------:|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| basis_file | string | Relative or absolute path to basis set file 

#### Optional ####
| Keyword | Type | Description | Default |
|------------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| bohr_coordinates | boolean | Set to true when coordinates specified in  $molecule section are in bohr units.  | False  |
| cart_bf | string | Controls the use of pure or Cartesian angular forms of GTOs. The letters corresponding to the angular momenta listed in this field will be expanded in cartesians, those not listed will be expanded in pure GTOs. For example, "pf" means p and f-type functions will be cartesian, and all others will be pure harmonic. Meanwhile, "spdf" means that s,p,d,and f-type functions will all be cartesian. Use caution and make sure that this field is specified correctly for your basis set, otherwise the dimensions of the %CAP matrix computed by OpenCAP and the TDMs computed by the electronic structure package won't match and the program will fail.  | "" (all functions are pure) |

## Projected_CAP ##
This section allows one to specify the parameters used to define the complex absorbing potential, the data read in from the electronic structure calculation, and the grid used for numerical integration.

There are two key pieces of data which must be read in from the electronic structure package: the set of one particle transition densities between each pair of states (each has a dimension of NxN where N is the number of basis functions), and the zeroth order Hamiltonian (which has a dimension of MxM, where M is the number of electronic states computed by the package). In some cases, the zeroth order Hamiltonian can be read directly from the output file of the electronic structure package, but in other cases it may be easier to just specify the zeroth order hamiltonian manually (due to its small dimensionality). For the latter option, the 'h0_file' keyword must be used, and the value must point to a properly formatted file. 

The ab initio basis set must be read in from a file formatted in the "Psi4" format. We suggest downloading a basis set from the [basis set exchange](https://www.basissetexchange.org/), and modifying it to suit your purposes. One crucial detail is that the ordering of the primitives in the basis set must match that of the electronic structure packaage. As a sanity check, OpenCAP will compute the overlap matrix and compare it to the overlap matrix computed by the electronic package (if available). If they match, you're good to go. If not, you should try adjusting the ordering of the primitives based on the indices of the overlap matrix elements OpenCAP reports as differing. If you're having trouble, a good place to check is the basis set library of your electronic structure package to see how the primitives are ordered there.

A properly formatted zeroth order hamiltonian file can have one of two forms. The first is the diagonal form (where the zeroth order hamiltonian just consists of total energies on the diagonal, as is the case for EOM-CCSD), the word "Diagonal" must appear on the first line, followed by the energies, each on their own line. The second is the full form, in which the entire zeroth order hamiltonian is specified. The word "Full" must appear on the first line, followed by the matrix. The matrix is read in row-major order. See [examples/heff.in](https://github.com/gayverjr/opencap/blob/master/example/heff.in) for an example.

We currently support two types of complex absorbing potentials; quadratic box-type potentials, and a quadratic "Voronoi isosurface" potential. Both types of %CAP have their own set of parameters which must be defined. Numerical integration is performed with the aid of the [numgrid](https://github.com/dftlibs/numgrid) package, which allocates the grid based on the molecular geometry and basis set. By default, a radial grid with precision 1e-14 is used, with an angular grid of 590 points. Both can optionally be adjusted to increase or decrease the size of the grid as needed.

#### Required ####
| Keyword | Type | Description | Valid options |
|:-------------:|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| cap_type | string | Type of absorbing potential| 'Box' or 'Voronoi'|
| Method | string | Name of electronic structure method | For OpenMolcas: MS-CASPT2, XMS-CASPT2, QD-NEVPT2 |
| Package | string | Name of electronic structure package | OpenMolcas |
| Nstates | int | Number of electronic states | 2 or greater |

#### OpenMolcas Interface ####
| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
| rassi_h5 | string | Relative or absolute path to rassi.h5 file containing TDMs.  |
| molcas_output | string | Relative or absolute path to OpenMolcas output file.  Required when h0_file keyword is unspecified.  |

#### Box CAP ####
| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
| CAP_X | float | Onset of %CAP in x-direction. Specify in bohr units regardless of whether 'bohr_coordinates' keyword is specified in $system section.|
| CAP_Y | float | Onset of %CAP in y-direction. Specify in bohr units regardless of whether 'bohr_coordinates' keyword is specified in $system section.|
| CAP_Z | float | Onset of %CAP in z-direction. Specify in bohr units regardless of whether 'bohr_coordinates' keyword is specified in $system section.|

#### Voronoi CAP ####
| Keyword | Type | Description |
|---------------|--------|------------------------------------------------------------------------------------------------------|
|r_cut | float | Cutoff radius for Voronoi %CAP. Specify in bohr units regardless of whether 'bohr_coordinates' keyword is specified in $system section.|


#### Optional ####
| Keyword | Type | Description | Default |
|------------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| radial_precision| int | Radial precision for numerical integration grid. A precision of 1x10^(-N), where N is the value specified is used.  | 14  |
| angular_points | int | Number of angular ponts used for the grid. See  https://github.com/dftlibs/numgrid for allowed numbers of points. | 590 |
| h0_file | string | Relative or absolute path to formatted zeroth order Hamiltonian file.  Required when output from electronic structure package (e.g. molcas_output) is unspecified. | "" |


