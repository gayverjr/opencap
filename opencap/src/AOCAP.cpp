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
}

// Smooth Voronoi CAP
// Thommas Sommerfeld and Masahiro Ehara
// DOI: 10.1021/acs.jctc.5b00465
void AOCAP::eval_voronoi_cap(double* x, double* y, double* z, double *grid_w, int num_points,Eigen::VectorXd &cap_values)
{
    #pragma omp parallel for
	for(size_t i=0;i<num_points;i++)
	{
		double atom_distances[num_atoms];
        double r_closest=1000;
		//find r closest and fill up our distances array
		for(size_t j=0;j<num_atoms;j++)
		{
			if(atoms[j].Z!=0)
			{
				  double dist_x= (x[i]-atoms[j].coords[0]) * (x[i]-atoms[j].coords[0]);
				  double dist_y= (y[i]-atoms[j].coords[1]) * (y[i]-atoms[j].coords[1]);
				  double dist_z= (z[i]-atoms[j].coords[2]) * (z[i]-atoms[j].coords[2]);
				  double dist=sqrt(dist_x+dist_y+dist_z);
				  if(dist<r_closest || j==0)
					 r_closest=dist;
				  atom_distances[j]=dist;
			}
		}
		double weights[num_atoms];
		for(size_t q=0;q<num_atoms;q++)
		{
		  if(atoms[q].Z!=0)
          {
              double weight = atom_distances[q]*atom_distances[q]-r_closest*r_closest+1;
              weights[q] = 1/(weight*weight);
          }
		  else
			  weights[q]=0;
		}
		double numerator=0.0;
		double denominator=0.0;
		for(size_t k=0;k<num_atoms;k++)
		{
		  numerator+=atom_distances[k]*atom_distances[k]*weights[k];
		  denominator+=weights[k];
		}
		double r=sqrt(numerator/denominator);
		double result=0;
		if(r>r_cut)
		  result = (r-r_cut)*(r-r_cut);
		cap_values(i) = result * grid_w[i];
	}
}

void AOCAP::eval_box_cap(double* x, double* y, double* z, double *grid_w, int num_points,Eigen::VectorXd &cap_values)
{
	for(size_t i=0;i<num_points;i++)
	{
		double result = 0;
		if(abs(x[i])>cap_x)
		 result += (abs(x[i])-cap_x) * (abs(x[i])-cap_x);
		if(abs(y[i])>cap_y)
		 result += (abs(y[i])-cap_y) * (abs(y[i])-cap_y);
		if(abs(z[i])>cap_z)
		 result += (abs(z[i])-cap_z) * (abs(z[i])-cap_z);
		cap_values(i) = result * grid_w[i];
	}
}

void AOCAP::eval_pot(double* x, double* y, double* z, double *grid_w, int num_points,Eigen::VectorXd &cap_values)
{
	if(compare_strings(cap_type, "box"))
		 eval_box_cap(x,y,z,grid_w,num_points,cap_values);
	else if (compare_strings(cap_type,"voronoi"))
		 eval_voronoi_cap(x,y,z,grid_w,num_points,cap_values);;
}

void AOCAP::compute_ao_cap_mat(Eigen::MatrixXd &cap_mat, BasisSet bs)
{
    std::cout << "Calculating CAP matrix in AO basis using " << std::to_string(omp_get_max_threads()) << " threads." << std::endl;
    std::cout << std::setprecision(2) << std::scientific  << "Radial precision: " << radial_precision
              << " Angular points: " << angular_points << std::endl;
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
    //double radial_precision = radial_precision;
    int min_num_angular_points = angular_points;
    int max_num_angular_points = angular_points;
	for(size_t i=0;i<num_atoms;i++)
	{
        
        // check radial parameters
        double alpha_max = bs.alpha_max(atoms[i]);
        std::vector<double> alpha_min = bs.alpha_min(atoms[i]);
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
                }
                else
                {
                    h = std::min(h,
                                 get_h(radial_precision, l, 0.1 * (r_outer - r_inner)));
                    if(r_outer < h)
                        opencap_throw("Error: r_outer < h, grid cannot be allocated for this basis.");
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
    	std::cout << "Grid for atom " + std::to_string(i+1) + " has: " << num_radial_points << " radial points, "
                  << angular_points << " angular points, " << num_points << " total points." << std::endl;
		evaluate_grid_on_atom(cap_mat,bs,grid_x_bohr,grid_y_bohr,grid_z_bohr,grid_w,num_points);
	}
}

void AOCAP::evaluate_grid_on_atom(Eigen::MatrixXd &cap_mat,BasisSet bs,double* grid_x_bohr,
		double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points)
{
	Eigen::VectorXd cap_values; Eigen::MatrixXd bf_values;
	cap_values = Eigen::VectorXd::Zero(num_points);
	bf_values = Eigen::MatrixXd::Zero(num_points,bs.num_carts());
	eval_pot(grid_x_bohr,grid_y_bohr,grid_z_bohr,grid_w,num_points,cap_values);
	size_t bf_idx = 0;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		std::vector<std::array<size_t,3>> order = opencap_carts_ordering(my_shell.l);
		for(size_t j=0;j<my_shell.num_carts();j++)
		{
			std::array<size_t,3> cart = order[j];
			my_shell.evaluate_on_grid(grid_x_bohr,grid_y_bohr,grid_z_bohr,num_points,cart[0],cart[1],cart[2],bf_values.col(bf_idx));
			bf_idx++;
		}
	}
	Eigen::MatrixXd bf_prime;
	bf_prime =  Eigen::MatrixXd::Zero(num_points,bs.num_carts());
	for(size_t i=0;i<bf_prime.cols();i++)
		bf_prime.col(i) = bf_values.col(i).array()*cap_values.array();
	cap_mat+=bf_prime.transpose()*bf_values;
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

