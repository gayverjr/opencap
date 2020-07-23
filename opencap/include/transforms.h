 /*! \file transforms.h
     \brief Functions for transforming between cartesian and solid harmonic GTOs, and for normalization
     of contracted GTOs.
 */
#include "BasisSet.h"
#include <Eigen/Dense>

/* Transforms matrix from cartesian to spherical harmonic representation.
 */
void cart2spherical(Eigen::MatrixXd &cart_ints, Eigen::MatrixXd &spherical_ints, BasisSet bs);
/* Sets all cartesian overlaps to 1.
 */
void uniform_cart_norm(Eigen::MatrixXd &my_mat, BasisSet bs);
