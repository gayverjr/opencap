/*
 * gto_ordering.h
 *
 *  Created on: Apr 27, 2020
 *      Author: JG
 */

#ifndef GTO_ORDERING_H_
#define GTO_ORDERING_H_


#include "BasisSet.h"
#include <armadillo>


void to_molcas_ordering(arma::mat &libcap_mat,BasisSet bs,std::vector<Atom> geometry);
std::vector<std::array<size_t,3>> molcas_carts_ordering(Shell shell);
std::vector<int> molcas_harmonic_ordering(Shell shell);
std::vector<std::vector<std::vector<Shell>>> molcas_reorder_basis_set(BasisSet bs, std::vector<Atom> geometry);
void to_molden_ordering(arma::mat &libcap_mat,BasisSet bs);
void to_qchem_ordering(arma::mat &libcap_mat,BasisSet bs);
void molden_to_qchem_ordering(arma::mat &libcap_mat,BasisSet bs);
std::vector<std::array<size_t,3>> molden_carts_ordering(Shell shell);
std::vector<int> molden_harmonic_ordering(Shell shell);



#endif /* GTO_ORDERING_H_ */
