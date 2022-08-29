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

/*
 * AOCAP.cpp
 */

#include "AOCAP.h"
#include "grid_radial.h"
#include "bragg.h"
#include <numgrid.h>
#include <omp.h>
#include <chrono>
#include <ctime>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <limits>
#include "BasisSet.h"
#include "gto_ordering.h"
#include "opencap_exception.h"
#include "utils.h"
#include <typeinfo>
#include "boxcap.h"
#include "cap_types.h"

AOCAP::AOCAP(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	verify_cap_parameters(params);
	double capx,capy,capz,rcut,radial,angpts,int_thresh;
	std::stringstream capxss(params["cap_x"]);
	std::stringstream capyss(params["cap_y"]);
	std::stringstream capzss(params["cap_z"]);
	std::stringstream rcutss(params["r_cut"]);
	std::stringstream radialss(params["radial_precision"]);
	std::stringstream angularss(params["angular_points"]);
	std::stringstream threshss(params["thresh"]);
	capxss >> capx;
	capyss >> capy;
	capzss >> capz;
	rcutss >> rcut;
	radialss >> radial;
	angularss >> angpts;
	threshss >> int_thresh;
	//now fill class members
	cap_type = params["cap_type"];
	cap_x = capx;
	cap_y= capy;
	cap_z = capz;
	r_cut = rcut;
	radial_precision = pow(10,-1.0*radial);
	angular_points = angpts;
	thresh = pow(10,-1.0*int_thresh);
	atoms = geometry;
    num_atoms = atoms.size();
	if(compare_strings(cap_type,"box"))
		cap_func = box_cap(cap_x,cap_y,cap_z);
	else if(compare_strings(cap_type,"voronoi"))
		cap_func = voronoi_cap(r_cut,atoms);
	else
		opencap_throw("Error: Only box and voronoi caps are implemented natively.");
	if(params.find("do_numerical")!=params.end())
	{
		std::string num_int = params["do_numerical"];
		if(compare_strings(num_int,"true"))
			do_numerical = true;
		else if(compare_strings(num_int,"false") && compare_strings(cap_type,"box"))
			do_numerical = false;
		else
			do_numerical = true; 
	}
	else if(compare_strings(cap_type,"box"))
		do_numerical = false;
	else
		do_numerical = true;
}

AOCAP::AOCAP(std::vector<Atom> geometry,std::map<std::string, std::string> params,const std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
std::vector<double> &, std::vector<double> &)> &custom_cap_func)
{
	verify_cap_parameters(params);
	double radial,angpts;
	std::stringstream radialss(params["radial_precision"]);
	std::stringstream angularss(params["angular_points"]);
	radialss >> radial;
	angularss >> angpts;
	//now fill class members
	cap_type = "custom";
	cap_x = 0.0;
	cap_y= 0.0;
	cap_z = 0.0;
	r_cut = 0.0;
	radial_precision = pow(10,-1.0*radial);
	angular_points = angpts;
	atoms = geometry;
    num_atoms = atoms.size();
	do_numerical = true;
	cap_func = custom_cap_func;
}


void AOCAP::integrate_cap_numerical(Eigen::MatrixXd &cap_mat, BasisSet bs)
{
	// funny business happens with GIL if we try to parallelize over atoms when custom python function is used
	if(cap_type=="custom")
		omp_set_num_threads(1);
    std::cout << "Calculating CAP matrix in AO basis using " << std::to_string(omp_get_max_threads()) << " threads." << std::endl;
    std::cout << std::setprecision(2) << std::scientific  << "Radial precision: " << radial_precision
              << " Angular points: " << angular_points << std::endl;
	size_t num_atoms = atoms.size();
    double x_coords_bohr[num_atoms];
	double y_coords_bohr[num_atoms];
	double z_coords_bohr[num_atoms];
	int nuc_charges[num_atoms];
	for(size_t i=0;i<num_atoms;i++)
	{
		x_coords_bohr[i]=atoms[i].coords[0];
		y_coords_bohr[i]=atoms[i].coords[1];
		z_coords_bohr[i]=atoms[i].coords[2];
		nuc_charges[i]=atoms[i].Z;
		if (atoms[i].Z==0)
			nuc_charges[i]=1; //choose bragg radius for H for ghost atoms
	}
    int min_num_angular_points = angular_points;
    int max_num_angular_points = angular_points;
	#pragma omp parallel for 
	for(size_t i=0;i<num_atoms;i++)
	{
        // check parameters
        double alpha_max = bs.alpha_max(atoms[i]);
        std::vector<double> alpha_min = bs.alpha_min(atoms[i],thresh);
        double r_inner = get_r_inner(radial_precision,
                                     alpha_max * 2.0); // factor 2.0 to match DIRAC
        double h = std::numeric_limits<float>::max();
        double r_outer = 0.0;
        for (int l = 0; l <= bs.max_L(); l++)
        {
            if (alpha_min[l] > 0.0)
            {
                r_outer =
                std::max(r_outer,
                         get_r_outer(radial_precision,
                                     alpha_min[l],
                                     l,
                                     4.0 * get_bragg_angstrom(nuc_charges[i])));
                if(r_outer < r_inner)
                {
                    opencap_throw("Error: r_outer < r_inner, grid cannot be allocated for this basis.");
                    //std::cout << "Setting alpha min[l] to 0.01" << std::endl; 
                    //alpha_min[l]=0.01;
                }
                else
                {
                    h = std::min(h,
                                 get_h(radial_precision, l, 0.1 * (r_outer - r_inner)));
                    if(r_outer < h)
                    {
                        opencap_throw("Error: r_outer < h, grid cannot be allocated for this basis.");
                    }
                }
            }
        }
        context_t *context = numgrid_new_atom_grid(radial_precision,
		                                 min_num_angular_points,
		                                 max_num_angular_points,
		                                 nuc_charges[i],
		                                 bs.alpha_max(atoms[i]),
		                                 bs.max_L(),
		                                 alpha_min.data());
		int num_points = numgrid_get_num_grid_points(context);
        double *grid_x_bohr = new double[num_points];
        double *grid_y_bohr = new double[num_points];
        double *grid_z_bohr = new double[num_points];
        double *grid_w = new double[num_points];
		double *cap_values = new double[num_points];
        numgrid_get_grid(  context,
                           num_atoms,
                           i,
                           x_coords_bohr,
                           y_coords_bohr,
                           z_coords_bohr,
                           nuc_charges,
                           grid_x_bohr,
                           grid_y_bohr,
                           grid_z_bohr,
                           grid_w);
        int num_radial_points = numgrid_get_num_radial_grid_points(context);
		compute_cap_on_grid(cap_mat,bs,grid_x_bohr,grid_y_bohr,grid_z_bohr,grid_w,num_points);
	}
}

void AOCAP::compute_cap_on_grid(Eigen::MatrixXd &cap_mat,BasisSet bs,double* x, double* y, double* z, 
double *grid_w, int num_points)
{
	Eigen::VectorXd cap_vals; Eigen::MatrixXd bf_values;
	bf_values = Eigen::MatrixXd::Zero(num_points,bs.num_carts());
    std::vector<double> x_vec(x,x+num_points);
    std::vector<double> y_vec(y,y+num_points);
    std::vector<double> z_vec(z,z+num_points);
    std::vector<double> w_vec(grid_w,grid_w+num_points);
    std::vector<double> cap_vec = cap_func(x_vec,y_vec,z_vec,w_vec);
	size_t bf_idx = 0;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		std::vector<std::array<size_t,3>> order = opencap_carts_ordering(my_shell.l);
		for(size_t j=0;j<my_shell.num_carts();j++)
		{
			std::array<size_t,3> cart = order[j];
			my_shell.evaluate_on_grid(x,y,z,num_points,cart[0],cart[1],cart[2],bf_values.col(bf_idx));
			bf_idx++;
		}
	}
    cap_vals = Eigen::Map<Eigen::VectorXd>(cap_vec.data(),cap_vec.size());
	Eigen::MatrixXd bf_prime;
	bf_prime =  Eigen::MatrixXd::Zero(num_points,bs.num_carts());
	for(size_t i=0;i<bf_prime.cols();i++)
		bf_prime.col(i) = bf_values.col(i).array()*cap_vals.array();
	cap_mat+=bf_prime.transpose()*bf_values;
}


void AOCAP::compute_ao_cap_mat(Eigen::MatrixXd &cap_mat, BasisSet &bs)
{
	if(!do_numerical)
	{
		std::cout << "CAP integrals will be computed analytically." << std::endl;
		eval_box_cap_analytical(cap_mat,bs);
		return;
	}
	else
		integrate_cap_numerical(cap_mat,bs);
}

void AOCAP::eval_box_cap_analytical(Eigen::MatrixXd &cap_mat, BasisSet &bs)
{
	double boxlength[3] = {cap_x,cap_y,cap_z};
	size_t bf1_idx = 0;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell shell1 = bs.basis[i];
		std::vector<std::array<size_t,3>> order1 = opencap_carts_ordering(shell1.l);
		for(size_t j=0;j<shell1.num_carts();j++)
		{
			std::array<size_t,3> l1 = order1[j];
			size_t bf2_idx = 0;
			for (size_t q=0;q<bs.basis.size();q++)
			{
				Shell shell2 = bs.basis[q];
				std::vector<std::array<size_t,3>> order2 = opencap_carts_ordering(shell2.l);
				for(size_t k=0;k<shell2.num_carts();k++)
				{
					std::array<size_t,3> l2 = order2[k];
					cap_mat(bf1_idx+j,bf2_idx+k) = integrate_box_cap(shell1,shell2,l1,l2,boxlength,thresh);
				}
				bf2_idx += shell2.num_carts();
			}
		}
		bf1_idx += shell1.num_carts();
	}
}


void AOCAP::verify_cap_parameters(std::map<std::string,std::string> &parameters)
{
	std::vector<std::string> missing_keys;
	if(parameters.find("cap_type")==parameters.end())
		opencap_throw("Error: Missing cap_type keyword.");
	if(compare_strings(parameters["cap_type"],"box"))
	{
		if(parameters.find("cap_x")==parameters.end())
			missing_keys.push_back("cap_x");
		if(parameters.find("cap_y")==parameters.end())
			missing_keys.push_back("cap_y");
		if (parameters.find("cap_z")==parameters.end())
			missing_keys.push_back("cap_z");
	}
	else if (compare_strings(parameters["cap_type"],"voronoi"))
	{
		if(parameters.find("r_cut")==parameters.end())
			missing_keys.push_back("r_cut");
	}
	else if(compare_strings(parameters["cap_type"],"custom"))
		;
	else
		opencap_throw("Error: only box and voronoi CAPs supported.");
	if(missing_keys.size()!=0)
	{
		std::string error_str = "Missing CAP keywords: ";
		for (auto key: missing_keys)
			error_str+=key+" ";
		opencap_throw(error_str);
	}
	std::map<std::string, std::string> defaults = {{"radial_precision", "16"}, {"angular_points", "590"}, {"thresh", "7"}};
	for (const auto &pair:defaults)
	{
		if(parameters.find(pair.first)==parameters.end())
			parameters[pair.first]=pair.second;
	}
}

