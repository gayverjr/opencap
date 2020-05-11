# OpenCAP input structure #

An OpenCAP input file has four sections; $geometry, $system, $cap_parameters, and $job. 
Each key/value pair must be specified by a newline, and each section must end with a 
$end command. For a simple example, refer to test.in the examples directory. We 
outline the valid keyword/value pairs for each of the sections below.

### Geometry ###
Molecular geometries are specified inline as follows:

atom x-coordinate y-coordinate z-coordinate

By default, coordinates are assumed to be in angstroms. To specify in bohr units, one must
set bohr_coordinates to true in the $system section.

### System ###
The system section specifies data to be read in from the electronic structure package, and data specifying the ab initio basis set.

We currently support an interface with the OpenMolcas quantum chemistry package, which specifies in multi-reference wave functions. 

#### Required ####

| Keyword | Type | Description | Valid options |
|:-------------:|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| Method | string | Name of electronic structure method | For OpenMolcas: MS-CASPT2, XMS-CASPT2, QD-NEVPT2 |
| Package | string | Name of electronic structure package | OpenMolcas |
| Nstates | int | Number of electronic states | 2 or greater |
| basis_file | string | Relative or absolute path to basis set file | N/A |
| rassi_h5 | string | Relative or absolute path to rassi.h5 file containing TDMs.  Only required for OpenMolcas interface. | N/A |
| molcas_output | string | Relative or absolute path to OpenMolcas output file.  Required when Package=OpenMolcas and h0_file keyword is unspecified.  | N/A |
| h0_file | string | Relative or absolute path to formatted zeroth order Hamiltonian file.  Required when output from electronic structure package (e.g. molcas_output) is unspecified. | N/A |
| basis_set | string | Title of basis set. Currently only "gen" is supported. | "gen" |

#### Optional ####

| Keyword | Type | Description | Default |
|------------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------|
| bohr_coordinates | boolean | Set to true when coordinates specified in  $molecule section are in bohr units.  | False  |
| cart_bf | string | Controls the use of pure or Cartesian angular forms of GTOs. The letters corresponding to the angular momenta listed in this field will be expanded in cartesians, those not listed will be expanded in pure GTOs. For example, "pf" means p and f-type functions will be cartesian, and all others will be pure harmonic. For comparison, "spdf" means that s,p,d,and f-type functions will all be cartesian. Use caution and make sure that this field is specified correctly for your basis set, otherwise the dimensions of the CAP matrix computed by OpenCAP and the TDMs computed by the electronic structure package won't match.  | "" (all functions are pure) |
|  |  |  |  |
