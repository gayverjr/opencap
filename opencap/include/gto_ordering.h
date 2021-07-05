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

/*! \file gto_ordering.h
     \brief Functions pertaining to ordering of Gaussian type orbitals (GTOs).
 */
#ifndef GTO_ORDERING_H_
#define GTO_ORDERING_H_

#include <Eigen/Dense>

#include "BasisSet.h"


/** Ordering of cartesian gaussians in molden ordering.
*/
std::vector<std::array<size_t,3>> opencap_carts_ordering(int l);
/** Ordering of harmonic gaussians in molden ordering.
 */
std::vector<int> opencap_harmonic_ordering(int l);
/* Get list of basis function IDs in OpenMolcas ordering.
 */
std::vector<bf_id> get_molcas_ids(BasisSet &bs,std::string rassi_filename);
/* Get list of basis function IDs in Q-Chem ordering.
 */
std::vector<bf_id> get_qchem_ids(BasisSet &bs);
/* Get list of basis function IDs in Psi4 ordering.
 */
std::vector<bf_id> get_psi4_ids(BasisSet &bs);
/* Get list of basis function IDs in Pyscf ordering.
 */
std::vector<bf_id> get_pyscf_ids(BasisSet &bs);
/* Re-orders matrix to opencap ordering of GTOs, given a set of input basis function ids
 * which describe the ordering of the original basis set.
 */
void to_opencap_ordering(Eigen::MatrixXd &mat,BasisSet &bs,std::vector<bf_id> input_ids);
/* Re-orders a matrix from ordering of original_ids to the ordering given by input_ids.
 */
void reorder_matrix(Eigen::MatrixXd &mat,std::vector<bf_id> &original_ids,std::vector<bf_id> &final_ids);


#endif /* GTO_ORDERING_H_ */
