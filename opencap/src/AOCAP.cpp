/*Copyright (c) 2020 James Gayvert

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
#include "BasisSet.h"
#include <numgrid.h>
#include "utils.h"
#include "gto_ordering.h"
#include "AOCAP.h"
#include <chrono>
#include <ctime>
#include <thread>
#include <vector>
#include <iostream>
#include "opencap_exception.h"
#include <Eigen/Dense>

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
}

// Smooth Voronoi CAP
// Thommas Sommerfeld and Masahiro Ehara
// DOI: 10.1021/acs.jctc.5b00465
double AOCAP::eval_voronoi_cap(double x, double y, double z)
{
    double atom_distances[atoms.size()];
    double r_closest=1000.0;
    //find r closest and fill up our distances array
    for(size_t j=0;j<atoms.size();j++)
    {
        if(atoms[j].Z!=0)
        {
			  double dist_x= (x-atoms[j].coords[0]) * (x-atoms[j].coords[0]);
			  double dist_y= (y-atoms[j].coords[1]) * (y-atoms[j].coords[1]);
			  double dist_z= (z-atoms[j].coords[2]) * (z-atoms[j].coords[2]);
			  double dist=sqrt(dist_x+dist_y+dist_z);
			  if(dist<r_closest)
				 r_closest=dist;
			  atom_distances[j]=dist;
        }
    }
    double weights[atoms.size()];
    for(size_t j=0;j<atoms.size();j++)
    {
  	  if(atoms[j].Z!=0)
  		  weights[j]=1/pow((pow(atom_distances[j],2.0)-pow(r_closest,2.0) +1),2.0);
  	  else
  		  weights[j]=0;
    }
    double numerator=0.0;
    double denominator=0.0;
    for(size_t j=0;j<atoms.size();j++)
    {
  	  numerator+=atom_distances[j]*atom_distances[j]*weights[j];
  	  denominator+=weights[j];
    }
    double r=sqrt(numerator/denominator);
    if(r<r_cut)
  	  return 0;
    else
  	  return (r-r_cut)*(r-r_cut);
}

double AOCAP::eval_box_cap(double x, double y, double z)
{
    double result = 0;
    if(abs(x)>cap_x)
   	 result += (abs(x)-cap_x) * (abs(x)-cap_x);
    if(abs(y)>cap_y)
   	 result += (abs(y)-cap_y) * (abs(y)-cap_y);
    if(abs(z)>cap_z)
   	 result += (abs(z)-cap_z) * (abs(z)-cap_z);
    return result;
}

double AOCAP::eval_pot(double x, double y, double z)
{
	if(compare_strings(cap_type, "box"))
		return eval_box_cap(x,y,z);
	else if (compare_strings(cap_type,"voronoi"))
		return eval_voronoi_cap(x,y,z);
	return 0;
}

void AOCAP::compute_ao_cap_mat(Eigen::MatrixXd &cap_mat, BasisSet bs)
{
	double x_coords_bohr[atoms.size()];
	double y_coords_bohr[atoms.size()];
	double z_coords_bohr[atoms.size()];
	int nuc_charges[atoms.size()];
	for(size_t i=0;i<atoms.size();i++)
	{
		x_coords_bohr[i]=atoms[i].coords[0];
		y_coords_bohr[i]=atoms[i].coords[1];
		z_coords_bohr[i]=atoms[i].coords[2];
		nuc_charges[i]=atoms[i].Z;
		if (atoms[i].Z==0)
			nuc_charges[i]=1; //choose bragg radius for H for ghost atoms
	}
    //double radial_precision = radial_precision;
    int min_num_angular_points = angular_points;
    int max_num_angular_points = angular_points;
	for(size_t i=0;i<atoms.size();i++)
	{
		//allocate and create grid
		context_t *context = numgrid_new_atom_grid(radial_precision,
		                                 min_num_angular_points,
		                                 max_num_angular_points,
		                                 nuc_charges[i],
		                                 bs.alpha_max(atoms[i]),
		                                 bs.max_L(),
		                                 &bs.alpha_min(atoms[i])[0]);
		int num_points = numgrid_get_num_grid_points(context);
        double *grid_x_bohr = new double[num_points];
        double *grid_y_bohr = new double[num_points];
        double *grid_z_bohr = new double[num_points];
        double *grid_w = new double[num_points];
        numgrid_get_grid(  context,
                           atoms.size(),
                           i,
                           x_coords_bohr,
                           y_coords_bohr,
                           z_coords_bohr,
                           nuc_charges,
                           grid_x_bohr,
                           grid_y_bohr,
                           grid_z_bohr,
                           grid_w);
		evaluate_grid_on_atom(cap_mat,bs,grid_x_bohr,grid_y_bohr,grid_z_bohr,grid_w,num_points);
	}

	//symmetrize
	for (size_t i=0;i<cap_mat.rows();i++)
	{
		for(size_t j=i+1;j<cap_mat.cols();j++)
			cap_mat(j,i) = cap_mat(i,j);
	}
}

void AOCAP::evaluate_grid_on_atom(Eigen::MatrixXd &cap_mat,BasisSet bs,double* grid_x_bohr,
		double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points)
{
	std::vector<std::vector<double>> bf_values;
	std::vector<double> cap_values (num_points);

	//pre-compute values of cap
	#pragma omp parallel for
	for (int i=0;i<num_points;i++)
		cap_values[i]= eval_pot(grid_x_bohr[i],grid_y_bohr[i],grid_z_bohr[i]);

	//pre-compute values of basis functions
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		std::vector<std::array<size_t,3>> order = opencap_carts_ordering(my_shell.l);
		for(size_t j=0;j<my_shell.num_carts();j++)
		{
			std::vector<double> vec(num_points);
			std::array<size_t,3> cart = order[j];
			#pragma omp parallel for
			for (int k=0;k<num_points;k++)
				vec[k]= my_shell.evaluate(grid_x_bohr[k],grid_y_bohr[k],grid_z_bohr[k],cart[0],cart[1],cart[2]);
			bf_values.push_back(vec);
		}
	}

    //Evaluating matrix elements
	#pragma omp parallel for
	for (size_t i=0;i<bs.num_carts();i++)
	{
		for(size_t j=i;j<bs.num_carts();j++)
		{
			for(int k=0;k<num_points;k++)
				cap_mat(i,j)+=grid_w[k]*cap_values[k]*bf_values[i][k]*bf_values[j][k];
		}
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

