/*
 * isodensity.h
 *
 *  Created on: Oct 16, 2020
 *      Author: JG
 */
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include "System.h"
#ifndef INCLUDE_ISODENSITY_H_
#define INCLUDE_ISODENSITY_H_


class isodensity{
public:

	Eigen::MatrixXd rho;
	System system;
	double isovalue;
	isodensity(System &my_sys, double iso, Eigen::MatrixXd my_rho);
	isodensity();
	double compute_distance_from_surface(std::vector<double> point);
	double eval_rho(std::vector<double> point);
	std::vector<double> numerical_gradient(std::vector<double> pt);
	std::vector<double> get_guess(std::vector<double> point);
};




#endif /* INCLUDE_ISODENSITY_H_ */
