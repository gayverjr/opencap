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
	void compute_cap_mat(arma::mat &Smat, BasisSet bs);

private:
	void num_overlap_block(Shell shell_a, Shell shell_b, arma::subview<double>&sub_mat,
			double* grid_x_bohr,double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points);
	double num_overlap_integral(Shell shell_a, std::array<size_t,3> a_cart, Shell shell_b,
			std::array<size_t,3> b_cart,double* grid_x_bohr,double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points);
	double eval_pot(double x, double y, double z);
	double eval_box_cap(double x, double y, double z);
	double eval_voronoi_cap(double x, double y, double z);
	size_t get_mat_idx(size_t bf_idx, BasisSet bs);
};



#endif /* CAP_H_ */
