/*
 * molcas_transform.cpp
 *
 *  Created on: Mar 15, 2020
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
#include "molcas_transform.h"


std::vector<std::array<size_t,3>> molcas_carts_ordering(Shell shell)
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
	else
		return {{0,0,0}};
}

std::vector<int> molcas_harmonic_ordering(Shell shell)
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

//molcas basis sets are specified a little differently, functions of the same subshell are grouped together
std::vector<std::vector<std::vector<Shell>>> molcas_reorder_basis_set(BasisSet bs, std::vector<Atom> geometry)
{
	//3D vector
	// Dimension 1: atom Shell belongs to
	// Dimension 2: subshell Shell belongs to
	// Dimension 3: actual Shell object
	std::vector<std::vector<std::vector<Shell>>> reordered_shells(geometry.size(),vector<std::vector<Shell>>(bs.max_L()+1));
	//group by atom and subshell
	for (auto shell:bs.basis)
	{
		//find which atom it belongs to
		size_t atm_idx = 0;
		while (shell.origin[2]!=geometry[atm_idx].coords[2])
			atm_idx++;
		reordered_shells[atm_idx][shell.l].push_back(shell);
	}
	return reordered_shells;
}

size_t find_matching_index(Shell shell, int angmom, BasisSet original_bs)
{
	std::vector<int> libcap_order = libcap_harmonic_ordering(shell);
	size_t matching_idx = 0;
	for (auto shell2:original_bs.basis)
	{
		if (shell==shell2)
			break;
		else
			matching_idx +=shell2.num_bf;
	}
	for (size_t i=0;i<libcap_order.size();i++)
	{
		if (libcap_order[i]==angmom)
			return matching_idx + i;
	}
	std::cout << "Something's gone wrong." << std::endl;
	return -1;
}


// matrix in libcap ordering --> matrix in molcas ordering
void to_molcas_ordering(arma::mat &libcap_mat, BasisSet bs, std::vector<Atom> geometry)
{
	std::vector<std::tuple<int,int>> swap_indices;
	std::vector<std::vector<std::vector<Shell>>> reordered_shells = molcas_reorder_basis_set(bs,geometry);
	size_t reordered_idx = 0;
	for (auto atm_group:reordered_shells)
	{
		for (auto subshell:atm_group)
		{
			if (!subshell.empty())
			{
				std::vector<int> molcas_order = molcas_harmonic_ordering(subshell[0]);
				for (auto angmom:molcas_order)
				{
					for(auto shell:subshell)
					{
						size_t original_idx = find_matching_index(shell,angmom,bs);
						swap_indices.push_back(std::make_tuple(original_idx, reordered_idx));
						reordered_idx++;
					}
				}
			}
		}
	}
	arma::mat per_mat(bs.Nbasis,bs.Nbasis);
	per_mat.zeros();
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	std::cout << std::endl;
	libcap_mat = per_mat.t()* libcap_mat * per_mat;
}







