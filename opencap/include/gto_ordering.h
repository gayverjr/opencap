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

/*! \file gto_ordering.h
     \brief Functions pertaining to ordering of Gaussian type orbitals (GTOs).
 */
#ifndef GTO_ORDERING_H_
#define GTO_ORDERING_H_
#include "BasisSet.h"
#include <Eigen/Dense>

/** Ordering of cartesian gaussians in opencap standard form.
*/
std::vector<std::array<size_t,3>> opencap_carts_ordering(int l);
/** Ordering of harmonic gaussians in opencap standard form.
 */
std::vector<int> opencap_harmonic_ordering(int l);
/* Get list of basis function IDs in OpenMolcas ordering.
 */
std::vector<bf_id> get_molcas_ids(BasisSet &bs,std::string rassi_filename);
/* Get list of basis function IDs in Molden ordering.
 */
std::vector<bf_id> get_molden_ids(BasisSet &bs);
/* Get list of basis function IDs in Q-Chem ordering.
 */
std::vector<bf_id> get_qchem_ids(BasisSet &bs);
/* Get list of basis function IDs in Pyscf ordering.
 */
std::vector<bf_id> get_pyscf_ids(BasisSet &bs);
/* Re-orders matrix to opencap ordering of GTOs, given a set of input basis function ids
 * which describe the ordering of the original basis set.
 */
void to_opencap_ordering(Eigen::MatrixXd &mat,BasisSet &bs,std::vector<bf_id> input_ids);


#endif /* GTO_ORDERING_H_ */
