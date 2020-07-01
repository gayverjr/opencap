 /*! \file qchem_interface.h
     \brief Functions pertaining to the Q-Chem interface.
 */
#ifndef QCHEM_INTERFACE_H_
#define QCHEM_INTERFACE_H_
#include <Eigen/Dense>
#include "BasisSet.h"
/*! Read in TDMs from Q-Chem .fchk file
    \param dmat_filename: file location of the Q-Chem .fchk file
    \param nstates: number of states
    \param num_bf: number of basis functions
    \return An array of size 2, each entry holding a 2D vector corresponding to the state densities.
     The first entry is the alpha densities, the second the beta densities.
*/
std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_in_dms_open_shell(std::string dmat_filename,size_t nstates,
		BasisSet bs);
std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_in_dms_closed_shell(std::string dmat_filename,size_t nstates,
		BasisSet bs);
Eigen::MatrixXd qchem_read_overlap(std::string dmat_filename, BasisSet bs);
/*! Read in energies from Q-Chem .fchk file
    \param nstates: number of states,
    \param method: electronic structure method method,
    \param output_file: qchem .fchk file location
    \return Eigen matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd read_qchem_energies(size_t nstates,std::string method,std::string output_file);

std::vector<Atom> read_atoms_from_fchk(std::string dmat_filename);
BasisSet read_basis_from_fchk(std::string dmat_filename, std::vector<Atom> atoms);
#endif /* QCHEM_INTERFACE_H_*/
