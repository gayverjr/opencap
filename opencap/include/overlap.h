 /*! \file overlap.h
     \brief Functions for computing the overlap matrix analytically.
 */
#ifndef OVERLAP_H_
#define OVERLAP_H_
#include <iostream>
#include <vector>
#include <math.h>
#include "Shell.h"
#include "BasisSet.h"
#include <Eigen/Dense>

/** Overlap between pair of cartesian gaussian functions
*/
double gauss_integral(double exp_a,std::array<size_t, 3> cart_a,std::array<double,3> coord_a,
		double exp_b,std::array<size_t, 3> cart_b, std::array<double,3> coord_b);
/** Overlap between pair of normalized, contracted gaussians
*/
double overlap_integral(Shell a, std::array<size_t,3> cart_a, Shell b, std::array<size_t,3> cart_);
/** Recursive formula for cartesian overlap distributions
*/
double mcmurchie_davidson(double exp_a, double exp_b, int t, int angmom_a, int angmom_b, double dist);
/** Computes group of integrals between pair of Shells
*/
Eigen::MatrixXd shell_overlap(Shell shell_a, Shell shell_b,Eigen::MatrixXd&sub_mat);
/** Computes overlap matrix
*/
void compute_analytical_overlap(BasisSet bs, Eigen::MatrixXd &Smat);

#endif /* OVERLAP_H_ */
