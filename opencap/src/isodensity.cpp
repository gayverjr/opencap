/*
 * isodensity.cpp
 *
 *  Created on: Oct 16, 2020
 *      Author: JG
 */
#include "isodensity.h"
#include "qchem_interface.h"
#include <pybind11/eigen.h>
#include "gto_ordering.h"
#include "BasisSet.h"
#include "transforms.h"
#include "Shell.h"
#include "utils.h"
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <iomanip>
#include <chrono>
#include <iostream>

isodensity::isodensity(System &my_sys, double iso, Eigen::MatrixXd my_rho)
{
	system = my_sys;
	isovalue = iso;
	rho = my_rho;
}

std::vector<double> isodensity::numerical_gradient(std::vector<double> pt)
{
	//central finite difference formula
	double h;
	h = 0.5;
	std::vector<double> forward,backward;
	//x-direction first
	h = pt[0]==0 ? sqrt(std::numeric_limits<double>::epsilon()) : sqrt(std::numeric_limits<double>::epsilon())*pt[0];
	forward = {pt[0]+h,pt[1],pt[2]};
	backward = {pt[0]-h,pt[1],pt[2]};
	double x_grad = (eval_rho(forward)-eval_rho(backward))/(2*h);
	h = pt[1]==0 ? sqrt(std::numeric_limits<double>::epsilon()) : sqrt(std::numeric_limits<double>::epsilon())*pt[1];
	forward = {pt[0],pt[1]+h,pt[2]};
	backward = {pt[0],pt[1]-h,pt[2]};
	double y_grad = (eval_rho(forward)-eval_rho(backward))/(2*h);
	h = pt[2]==0 ? sqrt(std::numeric_limits<double>::epsilon()) : sqrt(std::numeric_limits<double>::epsilon())*pt[2];
	forward = {pt[0],pt[1],pt[2]+h};
	backward = {pt[0],pt[1],pt[2]-h};
	double z_grad = (eval_rho(forward)-eval_rho(backward))/(2*h);
	return {x_grad,y_grad,z_grad};

}

double l2_norm(std::vector<double> const& u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

double dist(std::vector<double> const&u, std::vector<double> const &v)
{
	double dist = 0;
	for(size_t i=0;i<u.size();i++)
	{
		dist+= (u[i]-v[i])*(u[i]-v[i]);
	}
	return sqrt(dist);
}

std::vector<double> isodensity::get_guess(std::vector<double> point)
{
	double cur_val = eval_rho(point);
	if(cur_val > isovalue - isovalue*0.1)
		return point;


	std::vector<double> atm1 = {0.0,0.0,1.039};
	std::vector<double> atm2 = {0.0,0.0,-1.039};
	std::vector<double> COM = {0.0,0.0,0.0};

	std::vector<std::vector<double>> centers = {atm1,atm2,COM};
	std::vector<std::vector<double>> starting_points;
	std::vector<double> distances;

	for(auto &ctr:centers)
	{
		double x_dist = point[0] - ctr[0];
		double y_dist = point[1] - ctr[1];
		double z_dist = point[2] - ctr[2];
		std::vector<double> cur_point = point;
		double interval = 200;
		for(size_t i=1;i<=interval;i++)
		{
			cur_point[0] = cur_point[0]-x_dist/interval;
			cur_point[1] = cur_point[1]-y_dist/interval;
			cur_point[2] = cur_point[2]-z_dist/interval;
			cur_val = eval_rho(cur_point);
			if(cur_val>isovalue+isovalue*0.1)
			{
				//std::cout << "Overshot at:" << " x:" << cur_point[0] << " y:" << cur_point[1] << " z:" << cur_point[2] << " rho:" << cur_val << std::endl;
				//go back to previous point
				cur_point[0] = cur_point[0]+x_dist/interval;
				cur_point[1] = cur_point[1]+y_dist/interval;
				cur_point[2] = cur_point[2]+z_dist/interval;
				starting_points.push_back(cur_point);
				distances.push_back(dist(cur_point,point));
				break;
			}
		}
	}
	double min = 1E10;
	int min_idx = -1;
	for(size_t i=0;i<distances.size();i++)
	{
		if(distances[i]<min)
		{
			min=distances[i];
			min_idx = i;
		}
	}
	//std::cout << "Minimum was index:" << min_idx << std::endl;
	//should we ever get here? I don't think so
	return starting_points[min_idx];
}

double isodensity::compute_distance_from_surface(std::vector<double> point)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<double> cur_point = point;
	double cur_value = eval_rho(point);
	//std::cout << "Integration point:" << " x:" << cur_point[0] << " y:" << cur_point[1] << " z:" << cur_point[2] << " rho:" << cur_value << std::endl;
	size_t iter = 0;
    cur_point = get_guess(point);
	cur_value = eval_rho(cur_point);
	//std::cout << "Starting guess:" << " x:" << cur_point[0] << " y:" << cur_point[1] << " z:" << cur_point[2] << " rho:" << cur_value << std::endl;
	// newton's method
	while(cur_value < isovalue - isovalue*0.1)
	{
		std::vector<double> grad = numerical_gradient(cur_point);
		double grad_norm = l2_norm(grad);
		double disp = abs(cur_value - isovalue)/grad_norm;
		cur_point[0] = cur_point[0] + disp*grad[0];
		cur_point[1] = cur_point[1] + disp*grad[1];
		cur_point[2] = cur_point[2] + disp*grad[2];
		cur_value = eval_rho(cur_point);
		iter++;
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto total_time = std::chrono::duration<double>(stop-start).count();
	//std::cout << "Newton's method converged in " << total_time << " seconds and " << iter << " iterations." <<  std::endl;
	//std::cout << "Final point:" << " x:" << cur_point[0] << " y:" << cur_point[1] << " z:" << cur_point[2] << " rho:"<< cur_value <<  std::endl;
	double x_dist = point[0]-cur_point[0];
	double y_dist = point[1]-cur_point[1];
	double z_dist = point[2]-cur_point[2];
	double dist = sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist);
	//std::cout << "Calculated distance:" << dist << std::endl << std::endl;
	return dist;
}

Eigen::VectorXd eval_aos(std::vector<double> point,BasisSet bs)
{
	//
	Eigen::VectorXd cart_bf_values(bs.num_carts());
	Eigen::VectorXd spherical_bf_values(bs.Nbasis);
	spherical_bf_values = Eigen::VectorXd::Zero(bs.Nbasis);
	//calculate all values in cartesian basis and normalize appropriately
	size_t bf_idx = 0;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell shell = bs.basis[i];
		unsigned int cart_col_idx = 0;
		unsigned int sph_col_idx = 0;
		std::vector<std::array<size_t,3>> order = opencap_carts_ordering(shell.l);
		for(size_t j=0;j<shell.num_carts();j++)
		{
			std::array<size_t,3> cart = order[j];
			double scale = sqrt(fact2(2*shell.l-1)/fact2(2*cart[0]-1)
					/fact2(2*cart[1]-1)/fact2(2*cart[2]-1));
			cart_bf_values(bf_idx) = scale * shell.evaluate(point[0],point[1],point[2],
													cart[0],cart[1],cart[2]);
			bf_idx++;
		}
	}
	//rotate to spherical harmonic basis
	size_t cart_idx = 0;
	size_t sph_idx = 0;
	for(auto shell:bs.basis)
	{
		auto cart_seg = cart_bf_values.segment(cart_idx,shell.num_carts());
		spherical_bf_values.segment(sph_idx,shell.num_bf()) = get_trans_mat(shell)*cart_seg;
		cart_idx+=shell.num_carts();
		sph_idx+=shell.num_bf();
	}
	return spherical_bf_values;
}

double isodensity::eval_rho(std::vector<double> point)
{
	auto ao_vals = eval_aos(point,system.bs);
	return ao_vals.dot((rho * ao_vals).col(0));
}
