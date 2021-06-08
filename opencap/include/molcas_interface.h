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

/*! \file molcas_interface.h
     \brief Functions pertaining to the OpenMolcas interface.
 */

#include <Eigen/Dense>
#include "BasisSet.h"

#ifndef MOLCAS_INTERFACE_H_
#define MOLCAS_INTERFACE_H_

/*! Read in TDMs from OpenMolcas rassi.h5 file.
    \param alpha_opdms: MxM vector of vectors containing Eigen matrices for opdms
    \param beta_opdms: MxM vector of vectors containing Eigen matrices for opdms
    \param dmat_filename: file location of the rassi.h5 file.
    \param bs: BasisSet of the system
    \param nstates: Number of states
*/
void read_rassi_tdms(std::vector<std::vector<Eigen::MatrixXd>> &alpha_opdms,
		std::vector<std::vector<Eigen::MatrixXd>> &beta_opdms,
		std::string dmat_filename, BasisSet bs,size_t nstates);
/*! Read in overlap matrix from OpenMolcas rassi.h5 file.
    \param filename: String file location of the rassi.h5 file.
    \param bs: BasisSet of the system
    \return Eigen matrix containing the overlap matrix. Dimension is NxN, where N is number of basis functions.
*/
Eigen::MatrixXd read_rassi_overlap(std::string filename, BasisSet bs);
/*! Read in effective Hamiltonian from OpenMolcas output file.
    \param nstates: number of states
    \param filename: file location of the OpenMolcas output file.
    \return Eigen matrix containing the effective Hamiltonian matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd read_mscaspt2_heff(size_t nstates, std::string filename);
/*! Reads in basis set from rassi.h5 file.
 *  \param filename: file location of rassi.h5 file
 *  \param atoms: Ordered list of atoms in the system
 *  \return BasisSet specified by rassi.h5 file
*/
BasisSet read_basis_from_rassi(std::string filename,std::vector<Atom> atoms);
/*! Reads in geometry from rassi.h5 file.
 * \param filename: file location of rassi.h5 file
 * \return Geometry specified by rassi.h5 file
*/
std::vector<Atom> read_geometry_from_rassi(std::string filename);
/*! Read in effective Hamiltonian from OpenMolcas output file.
    \param nstates: number of states
    \param filename: file location of the OpenMolcas output file.
    \param method: pc-nevpt or sc-nevpt2
    \return Eigen matrix containing the effective Hamiltonian matrix. Dimension is MxM, where M is number of states.
*/
Eigen::MatrixXd read_nevpt2_heff(size_t nstates, std::string filename, std::string method);
#endif /* MOLCAS_INTERFACE_H_ */
