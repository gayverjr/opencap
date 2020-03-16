/*
 * molcas_transform.h
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */

#ifndef MOLCAS_TRANSFORM_H_
#define MOLCAS_TRANSFORM_H_

#include "BasisSet.h"
#include <armadillo>


void to_molcas_ordering(arma::mat &libcap_mat,BasisSet bs,std::vector<Atom> geometry);
std::vector<std::array<size_t,3>> molcas_carts_ordering(Shell shell);
std::vector<int> molcas_harmonic_ordering(Shell shell);
std::vector<Shell> molcas_reorder_basis_set(BasisSet bs, std::vector<Atom> geometry);



#endif /* MOLCAS_TRANSFORM_H_ */
