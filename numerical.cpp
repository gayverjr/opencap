/*
 * numerical.cpp
 *
 *  Created on: Feb 18, 2020
 *      Author: JG
 */
#include "BasisSet.h"
#include <armadillo>
#include <numgrid.h>
#include "utils.h"


double num_overlap_integral(Shell shell_a, std::array<size_t,3> a_cart, Shell shell_b,
		std::array<size_t,3> b_cart,double* grid_x_bohr,double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points)
{
	double total_integral = 0.0;
	#pragma omp parallel for reduction(+:total_integral)
	for(int i=0;i<num_points;i++)
	{
		double value=   grid_w[i] *
						shell_a.evaluate(grid_x_bohr[i],grid_y_bohr[i],grid_z_bohr[i],a_cart[0],a_cart[1],a_cart[2]) *
				        shell_b.evaluate(grid_x_bohr[i],grid_y_bohr[i],grid_z_bohr[i],b_cart[0],b_cart[1],b_cart[2]);
		total_integral+=value;
	}
	return total_integral;
}

void num_overlap_block(Shell shell_a, Shell shell_b, arma::subview<double>&sub_mat,
		double* grid_x_bohr,double *grid_y_bohr,double *grid_z_bohr,double *grid_w,int num_points)
{
	std::vector<std::array<size_t,3>> order_a = libcap_carts_ordering(shell_a);
	std::vector<std::array<size_t,3>> order_b = libcap_carts_ordering(shell_b);
	for(size_t i=0;i<shell_a.num_carts();i++)
	{
		std::array<size_t,3> a_cart = order_a[i];
		for(size_t j=0;j<shell_b.num_carts();j++)
		{
			std::array<size_t,3> b_cart = order_b[j];
			sub_mat(i,j) += num_overlap_integral(shell_a,a_cart, shell_b, b_cart,
					grid_x_bohr,grid_y_bohr,grid_z_bohr,grid_w,num_points);
		}
	}
}

void compute_numerical_overlap_mat(arma::mat &Smat, BasisSet bs,std::vector<Atom> atoms)
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
	}
    double radial_precision = 1.0e-10;
    int min_num_angular_points = 110;
    int max_num_angular_points = 110;
	for(size_t i=0;i<atoms.size();i++)
	{
		//allocate and create grid
		context_t *context = numgrid_new_atom_grid(radial_precision,
		                                 min_num_angular_points,
		                                 max_num_angular_points,
		                                 atoms[i].Z,
		                                 bs.alpha_max(atoms[i]),
		                                 bs.max_L(),
		                                 &bs.alpha_min(atoms[i])[0]);
		int num_points = numgrid_get_num_grid_points(context);
        double *grid_x_bohr = new double[num_points];
        double *grid_y_bohr = new double[num_points];
        double *grid_z_bohr = new double[num_points];
        double *grid_w = new double[num_points];
        std::cout << "Getting the grid for atom:" << i+1 << std::endl;
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
		unsigned int row_idx = 0;
		std::cout << "Finished allocating the grid for atom:" << i+1 << std::endl;
		std::cout << "Grid has " << num_points << " points." << std::endl;
		for(size_t i=0;i<bs.basis.size();i++)
		{
			Shell shell1 = bs.basis[i];
			unsigned int col_idx = 0;
			for(size_t j=0;j<bs.basis.size();j++)
			{
				Shell shell2= bs.basis[j];
				//view to block of the matrix corresponding to these two pairs of basis functions
				auto sub_mat = Smat.submat(row_idx,col_idx,
						row_idx+shell1.num_carts()-1,col_idx+shell2.num_carts()-1);
				num_overlap_block(shell1,shell2,sub_mat,grid_x_bohr,grid_y_bohr,grid_z_bohr,grid_w,num_points);
				col_idx += shell2.num_carts();
			}
			row_idx += shell1.num_carts();
		}
        delete[] grid_x_bohr;
        delete[] grid_y_bohr;
        delete[] grid_z_bohr;
        delete[] grid_w;
        numgrid_free_atom_grid(context);
	}
}

