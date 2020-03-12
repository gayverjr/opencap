/*
 * molden_transform.h
 *
 *  Created on: Mar 12, 2020
 *      Author: JG
 */

#ifndef MOLDEN_TRANSFORM_H_
#define MOLDEN_TRANSFORM_H_
#include "BasisSet.h"
#include <armadillo>


void to_molden_ordering(arma::mat &libcap_mat,BasisSet bs);
std::vector<std::array<size_t,3>> molden_carts_ordering(Shell shell);
std::vector<int> molden_harmonic_ordering(Shell shell);



#endif /* MOLDEN_TRANSFORM_H_ */
