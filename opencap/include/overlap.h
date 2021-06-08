/*Copyright (c) 2021 James Gayvert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*! \file overlap.h
     \brief Functions for computing the overlap matrix analytically.
 */
#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "BasisSet.h"


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
