 /*! \file CAP.h
     \brief Class for numerically integrating the %CAP matrix in AO basis.
 */
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include "Atom.h"
#include <armadillo>
#include "BasisSet.h"
#pragma once

/*! \brief Class for numerically integrating the %CAP matrix in AO basis.
 *
 */

class CAP
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
	CAP(std::vector<Atom> geometry,std::map<std::string, std::string> params);
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
	/** Geometry of molecular system.
	 */
	std::vector<Atom> atoms;
	/** Computes %CAP matrix in AO basis via numerical integration.
	 */
	void compute_cap_mat(arma::mat &cap_mat, BasisSet bs);

private:
	/** Evaluate potential at grid point.
	 */
	double eval_pot(double x, double y, double z);
	/** Evaluate box %CAP at grid point.
	 */
	double eval_box_cap(double x, double y, double z);
	/** Evaluate Voronoi %CAP at grid point.
	 */
	double eval_voronoi_cap(double x, double y, double z);
	/** Evaluate all points on grid for a given atom.
	 */
	void evaluate_grid_on_atom(arma::mat &cap_mat,BasisSet bs,double* grid_x_bohr,
			double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points);
};
