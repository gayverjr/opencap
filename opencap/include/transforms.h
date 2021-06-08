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

/*! \file transforms.h
     \brief Functions for transforming between cartesian and solid harmonic GTOs, and for normalization
     of contracted GTOs.
 */
#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

#include <Eigen/Dense>


/* Transforms matrix from cartesian to spherical harmonic representation.
 */
void cart2spherical(Eigen::MatrixXd &cart_ints, Eigen::MatrixXd &spherical_ints, BasisSet &bs);
/* Sets all cartesian overlaps to 1.
 */
void uniform_cart_norm(Eigen::MatrixXd &my_mat, BasisSet &bs);

#endif /* TRANSFORMS_H_ */
