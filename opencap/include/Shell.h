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
/*! \file Shell.h
     \brief Class for a group of atomic orbitals with the same angular momentum and
     exponents.
 */

#include <array>
#include <vector>

#pragma once
#include "Atom.h"
#include <Eigen/Dense>
using namespace std;

/*! \brief Class for a group of atomic orbitals with the same angular momentum and
     exponents.
 *
 *  A %shell is a group of atomic orbitals in the ab initio basis set which share the same orbital
 *  angular momentum quantum number (l)and are defined by the same group of gaussian primitves.
 *  A shell is defined by a set of exponents and expansion coefficients.  To reduce the number of
 *  basis functions and  reduce linear dependence, GTOs can be expressed as pure solid harmonics
 *  as opposed to cartesian GTOs.
 *
 */

class Shell {
public:
	/** Orbital angular momentum quantum number
	 */
	int l;
	/** True = pure solid harmonic GTOs, False = cartesian GTOs
	 */
	bool pure;
	/** List of exponents
	 */
	std::vector<double> exps;
	/** Expansion coefficients for contraction
	 */
	std::vector<double> coeffs;
	/** Number of primitive gaussians in the %Shell
	 */
	size_t num_prims;
	/** Center of gaussian orbitals
	 */
	std::array<double,3> origin;
	/** Number of cartesian basis functions
	 */
	size_t num_carts();
	/** Constructor from data read in from file
	 */
	Shell(int angmom, bool cart_flag,std::vector<double> exponents,
			std::vector<double> coefficients);
	/** Copy constructor with coordinates
	 */
	Shell(Shell s, std::array<double,3>center);
	/** Default constructor
	 */
	Shell();
	Shell(int angmom,std::array<double,3>center);
	/** Updates center of %Shell
	 */
	void update_coords(std::array<double,3> center);
    /** Evaluate basis function (with a given lx ly lz) at point x,y,z
     */
	double evaluate(double x, double y, double z,size_t lx,size_t ly, size_t lz);
    /** Evaluate basis function (with a given lx ly lz) on grid
     */
	void evaluate_on_grid(double* x, double* y, double* z, int num_points,size_t lx,size_t ly, size_t lz,Eigen::Ref<Eigen::VectorXd> v);
	/** Returns smallest exponent in %Shell
	 */
	std::vector<double> alpha_min();
	/** Returns largest exponent in %Shell
	 */
	double alpha_max();
	/** Checks if two shells are equivalent (same exps,coeffs,center)
	 */
    bool operator==(const Shell& other);
    void add_primitive(double exp,double coeff);
	/** Normalizes contraction coefficients so that self overlap = 1
	 */
	void normalize();
	/** Returns number of basis functions
	 */
	size_t num_bf();

private:


};
