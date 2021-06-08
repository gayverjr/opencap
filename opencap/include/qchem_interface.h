/*Copyright (c) 2021 James Gayvert

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
SOFTWARE.
*/

/*! \file qchem_interface.h
     \brief Functions pertaining to the Q-Chem interface.
 */
#ifndef QCHEM_INTERFACE_H_
#define QCHEM_INTERFACE_H_

#include <Eigen/Dense>
#include "BasisSet.h"

/*! Reads in overlap matrix from Q-Chem fchk file.
    \param dmat_filename: file location of the Q-Chem .fchk file
    \param bs: BasisSet of system
    \return Eigen matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd qchem_read_overlap(std::string dmat_filename, BasisSet bs);
/*! Read in energies from Q-Chem output file for EOM-EE/EA/IP methods
    \param nstates: number of states,
    \param output_file: qchem output file location
    \return Eigen matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd read_qchem_eom_energies(size_t nstates,std::string output_file);
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
/*! Read in TDMs from Q-Chem fchk file.
    \param alpha_opdms: MxM vector of vectors containing Eigen matrices for opdms
    \param beta_opdms: MxM vector of vectors containing Eigen matrices for opdms
    \param fchk_filename: file location of the rassi.h5 file.
    \param bs: BasisSet of the system
    \param nstates: Number of states
*/
void qchem_read_dms(std::vector<std::vector<Eigen::MatrixXd>> &alpha_dms,
		std::vector<std::vector<Eigen::MatrixXd>> &beta_dms,
		std::string fchk_filename,BasisSet &bs, size_t num_states);
/*! Read in energies from Q-Chem output file for EOM-EE/EA/IP methods
    \param nstates: number of states,
    \param output_file: qchem output file location
    \return Eigen matrix containing the diagonal H0 matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd read_qchem_tddft_energies(size_t nstates,std::string output_file);

#endif /* QCHEM_INTERFACE_H_*/
