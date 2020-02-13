/*
 * transforms.h
 *
 *  Created on: Feb 12, 2020
 *      Author: JG
 */
#include <armadillo>
#include "BasisSet.h"

void cart2spherical(arma::mat &cart_ints, BasisSet bs);
void uniform_cart_norm(arma::mat &my_mat, BasisSet bs);
double get_coeff(int L, int m, int lx, int ly, int lz);
