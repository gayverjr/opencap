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
#include "boxcap.h"
#include "numerical_integration.h"
#include "cap_types.h"

AOCAP::AOCAP(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	verify_cap_parameters(params);
	double capx,capy,capz,rcut,radial,angpts;
	std::stringstream capxss(params["cap_x"]);
	std::stringstream capyss(params["cap_y"]);
	std::stringstream capzss(params["cap_z"]);
	std::stringstream rcutss(params["r_cut"]);
	std::stringstream radialss(params["radial_precision"]);
	std::stringstream angularss(params["angular_points"]);
	capxss >> capx;
	capyss >> capy;
	capzss >> capz;
	rcutss >> rcut;
	radialss >> radial;
	angularss >> angpts;
	//now fill class members
	cap_type = params["cap_type"];
	cap_x = capx;
	cap_y= capy;
	cap_z = capz;
	r_cut = rcut;
	radial_precision = pow(10,-1.0*radial);
	angular_points = angpts;
	atoms = geometry;
    num_atoms = atoms.size();
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

void AOCAP::compute_ao_cap_mat(Eigen::MatrixXd &cap_mat, BasisSet &bs)
{
	if(!do_numerical)
	{
		std::cout << "CAP integrals will be computed analytically." << std::endl;
		eval_box_cap_analytical(cap_mat,bs);
		return;
	}
	else
	{
		std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
		std::vector<double> &, std::vector<double> &, int)> cap_func;
		if(compare_strings(cap_type,"box"))
			cap_func = box_cap(cap_x,cap_y,cap_z);
		else if(compare_strings(cap_type,"voronoi"))
			cap_func = voronoi_cap(r_cut,atoms);
		else
			opencap_throw("Error: Only box and voronoi caps are implemented natively.");
		integrate_cap_numerical(cap_mat,bs,atoms,radial_precision,angular_points,cap_func);
	}
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
					cap_mat(bf1_idx+j,bf2_idx+k) = integrate_box_cap(shell1,shell2,l1,l2,boxlength);
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
	else
		opencap_throw("Error: only box and voronoi CAPs supported.");
	if(missing_keys.size()!=0)
	{
		std::string error_str = "Missing CAP keywords: ";
		for (auto key: missing_keys)
			error_str+=key+" ";
		opencap_throw(error_str);
	}
	std::map<std::string, std::string> defaults = {{"radial_precision", "14"}, {"angular_points", "590"}};
	for (const auto &pair:defaults)
	{
		if(parameters.find(pair.first)==parameters.end())
			parameters[pair.first]=pair.second;
	}
}

