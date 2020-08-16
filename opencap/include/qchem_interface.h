/*Copyright (c) 2020 James Gayvert

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
std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_dms(std::string fchk_filename,BasisSet &bs);
#endif /* QCHEM_INTERFACE_H_*/
