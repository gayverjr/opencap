 /*! \file gto_ordering.h
     \brief Functions pertaining to ordering of Gaussian type orbitals (GTOs).
 */
#ifndef GTO_ORDERING_H_
#define GTO_ORDERING_H_
#include "BasisSet.h"
#include <armadillo>

/** Ordering of cartesian gaussians in opencap.
 */
std::vector<std::array<size_t,3>> opencap_carts_ordering(Shell shell);
/** Ordering of harmonic gaussians in opencap.
 */
std::vector<int> opencap_harmonic_ordering(Shell shell);
/** Re-order matrix from opencap ordering to molcas ordering.
 */
void to_molcas_ordering(arma::mat &opencap_mat,BasisSet bs,std::vector<Atom> geometry);
/** Ordering of cartesian gaussians in molcas.
 */
std::vector<std::array<size_t,3>> molcas_carts_ordering(Shell shell);
/** Ordering of harmonic gaussians in molcas.
 */
std::vector<int> molcas_harmonic_ordering(Shell shell);
/** Re-order shells in molcas format.
 */
std::vector<std::vector<std::vector<Shell>>> molcas_reorder_basis_set(BasisSet bs, std::vector<Atom> geometry);
/** Re-order matrix from opencap ordering to molden.
 */
void to_molden_ordering(arma::mat &opencap_mat,BasisSet bs);
/** Re-order matrix from opencap ordering to q-chem internal ordering.
 */
void to_qchem_ordering(arma::mat &opencap_mat,BasisSet bs);
/** Re-order matrix from molden ordering to q-chem internal ordering.
 */
void molden_to_qchem_ordering(arma::mat &opencap_mat,BasisSet bs);
/** Ordering of cartesian gaussians in molden format.
 */
std::vector<std::array<size_t,3>> molden_carts_ordering(Shell shell);
/** Ordering of harmonic gaussians in molden format.
 */
std::vector<int> molden_harmonic_ordering(Shell shell);



#endif /* GTO_ORDERING_H_ */
