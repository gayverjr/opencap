/*
 * Numerical_CAP.h
 *
 *  Created on: Apr 2, 2020
 *      Author: JG
 */
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include "Atom.h"
#include <armadillo>
#include "BasisSet.h"
#ifndef CAP_H_
#define CAP_H_
class CAP
{
public:
	double radial_precision;
	size_t angular_points;
	CAP(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	std::string cap_type;
	double cap_x;
	double cap_y;
	double cap_z;
	double r_cut;
	std::vector<Atom> atoms;
	void compute_cap_mat(arma::mat &cap_mat, BasisSet bs);

private:
	double eval_pot(double x, double y, double z);
	double eval_box_cap(double x, double y, double z);
	double eval_voronoi_cap(double x, double y, double z);
	void evaluate_grid_on_atom(arma::mat &cap_mat,BasisSet bs,double* grid_x_bohr,
			double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points);
};



#endif /* CAP_H_ */
