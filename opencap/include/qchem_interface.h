 /*! \file qchem_interface.h
     \brief Functions pertaining to the Q-Chem interface.
 */
#ifndef QCHEM_INTERFACE_H_
#define QCHEM_INTERFACE_H_
#include <Eigen/Dense>
#include "BasisSet.h"
/*! Read in TDMs from Q-Chem .fchk file for open-shell system (separate alpha and beta densities)
    \param dmat_filename: file location of the Q-Chem .fchk file
    \param nstates: number of states
    \param bs: BasisSet of system
    \return An array of size 2, each entry holding a 2D vector corresponding to the state densities.
     The first entry is the alpha densities, the second the beta densities.
*/
std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_in_dms_open_shell(std::string dmat_filename,size_t nstates,
		BasisSet bs);
/*! Read in TDMs from Q-Chem .fchk file for open-shell system (no distinction between alpha/beta densities)
    \param dmat_filename: file location of the Q-Chem .fchk file
    \param nstates: number of states
    \param bs: BasisSet of system
    \return An array of size 2, each entry holding a 2D vector corresponding to the state densities.
     The first entry is the alpha densities, the second the beta densities.
*/
std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_in_dms_closed_shell(std::string dmat_filename,size_t nstates,
		BasisSet bs);
/*! Reads in overlap matrix from Q-Chem fchk file.
    \param dmat_filename: file location of the Q-Chem .fchk file
    \param bs: BasisSet of system
    \return Eigen matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd qchem_read_overlap(std::string dmat_filename, BasisSet bs);
/*! Read in energies from Q-Chem output file for EOM-EE/EA/IP methods
    \param nstates: number of states,
    \param method: electronic structure method method,
    \param output_file: qchem .fchk file location
    \return Eigen matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd read_qchem_energies(size_t nstates,std::string method,std::string output_file);
/* Reads geometry from Q-Chem fchk file.
  \param fchk_filename: qchem .fchk file location
  \return atoms: List of atoms in system
 */
std::vector<Atom> read_geometry_from_fchk(std::string fchk_filename);
/* Reads basis set from Q-Chem fchk file.
     \param fchk_filename: qchem .fchk file location
     \return bs: BasisSet of system
 */
BasisSet read_basis_from_fchk(std::string fchk_filename, std::vector<Atom> atoms);
#endif /* QCHEM_INTERFACE_H_*/
