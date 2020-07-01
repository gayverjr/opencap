 /*! \file transforms.h
     \brief Functions for transforming between cartesian and solid harmonic GTOs, and for normalization
     of contracted GTOs.
 */
#include "BasisSet.h"
#include <Eigen/Dense>

void cart2spherical(Eigen::MatrixXd &cart_ints, Eigen::MatrixXd &spherical_ints, BasisSet bs);
void uniform_cart_norm(Eigen::MatrixXd &my_mat, BasisSet bs);
double get_coeff(int L, int m, int lx, int ly, int lz);
Eigen::MatrixXd transform_block(Shell shell1, Shell shell2, Eigen::MatrixXd cart_block);
Eigen::MatrixXd get_trans_mat(Shell shell,std::string ordering);
