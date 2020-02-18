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
void transform_block(Shell shell1, Shell shell2, arma::subview<double>&cart_sub_mat,arma::subview<double>&sph_sub_mat);
arma::mat get_trans_mat(Shell shell);
