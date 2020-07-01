 /*! \file gto_ordering.h
     \brief Functions pertaining to ordering of Gaussian type orbitals (GTOs).
 */
#ifndef GTO_ORDERING_H_
#define GTO_ORDERING_H_
#include "BasisSet.h"
#include <Eigen/Dense>

/** Ordering of cartesian gaussians in opencap.
 */
std::vector<std::array<size_t,3>> opencap_carts_ordering(Shell shell);
/** Ordering of harmonic gaussians in opencap.
 */
std::vector<int> opencap_harmonic_ordering(Shell shell);
/** Re-order matrix from opencap ordering to molcas ordering.
 */
void to_molcas_ordering(Eigen::MatrixXd &opencap_mat,BasisSet bs,std::vector<Atom> geometry);
/** Re-order matrix from opencap ordering to pyscf ordering.
 */
void to_pyscf_ordering(Eigen::MatrixXd &opencap_mat,BasisSet bs);
/** Ordering of cartesian gaussians in molcas.
 */
std::vector<std::array<size_t,3>> molcas_carts_ordering(Shell shell);
/** Ordering of harmonic gaussians in molcas.
 */
std::vector<int> molcas_harmonic_ordering(Shell shell);
/** Ordering of cartesian gaussians in molden format.
 */
std::vector<std::array<size_t,3>> molden_carts_ordering(Shell shell);
/** Ordering of harmonic gaussians in molden format.
 */
std::vector<int> molden_harmonic_ordering(Shell shell);
/** Ordering of harmonic gaussians in pyscf format.
 */
std::vector<int> pyscf_harmonic_ordering(Shell shell);
/** Ordering of cartesian gaussians in pyscf format.
 */
std::vector<bf_id> get_molcas_ids(BasisSet bs,std::string rassi_filename);
std::vector<bf_id> get_molden_ids(BasisSet bs);
std::vector<bf_id> get_pyscf_ids(BasisSet bs);
void to_opencap_ordering(Eigen::MatrixXd &mat,BasisSet bs,std::vector<bf_id> input_ids);



#endif /* GTO_ORDERING_H_ */
