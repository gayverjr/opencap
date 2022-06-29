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

/*! \file AOCAP.h
     \brief Class for numerically integrating the %CAP matrix in AO basis.
 */

#include <cstdlib>
#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

#include "Atom.h"
#include "BasisSet.h"

#pragma once

/*! \brief Class for numerically integrating the %CAP matrix in AO basis.
 *
 */

class AOCAP
{
public:
	/** Radial precision of numerical grid. Default is 1.0e-14
	 */
	double radial_precision;
	/** Number of angular points on grid. Default is 590. See https://github.com/dftlibs/numgrid
	 *  for allowed number of points.
	 */
	size_t angular_points;
	/** Constructs %CAP object from geometry and CAP parameters.
	 */
	AOCAP(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	AOCAP(std::vector<Atom> geometry,std::map<std::string, std::string> params,const std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
std::vector<double> &, std::vector<double> &)> &cap_func);
	/** Default construct, does nothing
	*/
	AOCAP(){cap_x=0.0;cap_y=0.0;cap_z=0.0;r_cut=0.0;do_numerical=true;num_atoms=0;};
	/** Type of %CAP. Can be Voronoi or Box CAP.
	 */
	std::string cap_type;
	/** Onset of box %CAP in X direction. Specify in bohr units.
	 */
	double cap_x;
	/** Onset of box %CAP in Y direction. Specify in bohr units.
	 */
	double cap_y;
	/** Onset of box %CAP in Z direction. Specify in bohr units.
	 */
	double cap_z;
	/** Cutoff radius of Voronoi %CAP. Specify in bohr units.
	 */
	double r_cut;
	/** Threshold for exponents of GTOs. Exponents smaller than thresh are excluded from integration. Intended to alleviate numerical issues associated with fake IP orbital.
	*/
	double thresh;
	/** Geometry of molecular system.
	 */
	std::vector<Atom> atoms;
	/** Computes %CAP matrix in AO basis via numerical integration.
	 */
	void compute_ao_cap_mat(Eigen::MatrixXd &cap_mat, BasisSet &bs);
    /** Number of atoms
     */
    size_t num_atoms;
	bool do_numerical;
	std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
		std::vector<double> &, std::vector<double> &)> cap_func;
	void compute_cap_on_grid(Eigen::MatrixXd &cap_mat,BasisSet bs,double* x, double* y, double* z, double *grid_w, int num_points);

private:
	/** Checks whether specified CAP is valid.
	 */
	void verify_cap_parameters(std::map<std::string,std::string> &parameters);
	void eval_box_cap_analytical(Eigen::MatrixXd &cap_mat, BasisSet &bs);
	void integrate_cap_numerical(Eigen::MatrixXd &cap_mat, BasisSet bs);
};
