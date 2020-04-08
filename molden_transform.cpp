/*
 * molden_transform.cpp
 *
 *  Created on: Mar 10, 2020
 *      Author: JG
 */
#include <armadillo>
#include "BasisSet.h"
#include "Shell.h"
#include "utils.h"
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include <tuple>
#include "molden_transform.h"


std::vector<std::array<size_t,3>> molden_carts_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {{0,0,0}};
	//p
	if(shell.l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (shell.l==2)
		return {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
	//f
	else if (shell.l==3)
		return {{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},{2,0,1},{1,0,2},
				{0,1,2},{0,2,1},{1,1,1}};
	//g
	else if(shell.l==4)
		return {{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},{1,3,0},{0,3,1},
		        {1,0,3},{0,1,3},{2,2,0},{2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,0,2}};
	else
		return {{0,0,0}};
}

std::vector<int> molden_harmonic_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {0};
	//p
	if(shell.l==1)
		return {1,-1,0};
	//d
	else if (shell.l==2)
		return {0,1,-1,2,-2};
	//f
	else if (shell.l==3)
		return {0,1,-1,2,-2,3,-3};
	//g
	else if (shell.l==4)
		return {0,1,-1,2,-2,3,-3,4,-4};
	else
		return {0,0,0};
}

std::vector<int> qchem_harmonic_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {0};
	//p
	if(shell.l==1)
		return {1,-1,0};
	//d
	else if (shell.l==2)
		return {-2,-1,0,1,2};
	//f
	else if (shell.l==3)
		return {-3,-2,-1,0,1,2,3};
	//g
	else if (shell.l==4)
		return {-4,-3,-2,-1,0,1,2,3,4};
	else
		return {0,0,0};
}

// matrix in libcap ordering --> matrix in molden ordering
void to_molden_ordering(arma::mat &libcap_mat, BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> libcap_order = libcap_harmonic_ordering(shell);
			std::vector<int> molden_order = molden_harmonic_ordering(shell);
			for (size_t i=0;i<libcap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(libcap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> libcap_order = libcap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(shell);
			for (size_t i=0;i<libcap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(libcap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	arma::mat per_mat(bs.Nbasis,bs.Nbasis);
	per_mat.zeros();
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	libcap_mat = per_mat.t()* libcap_mat * per_mat;
}

// matrix in libcap ordering --> matrix in qchem ordering
void to_qchem_ordering(arma::mat &libcap_mat, BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> libcap_order = libcap_harmonic_ordering(shell);
			std::vector<int> molden_order = qchem_harmonic_ordering(shell);
			for (size_t i=0;i<libcap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(libcap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> libcap_order = libcap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(shell);
			for (size_t i=0;i<libcap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(libcap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	arma::mat per_mat(bs.Nbasis,bs.Nbasis);
	per_mat.zeros();
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	libcap_mat = per_mat.t()* libcap_mat * per_mat;
}

// matrix in libcap ordering --> matrix in qchem ordering
void molden_to_qchem_ordering(arma::mat &libcap_mat, BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> libcap_order = molden_harmonic_ordering(shell);
			std::vector<int> molden_order = qchem_harmonic_ordering(shell);
			for (size_t i=0;i<libcap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(libcap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> libcap_order = libcap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(shell);
			for (size_t i=0;i<libcap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(libcap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	arma::mat per_mat(bs.Nbasis,bs.Nbasis);
	per_mat.zeros();
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	libcap_mat = per_mat.t()* libcap_mat * per_mat;
}
