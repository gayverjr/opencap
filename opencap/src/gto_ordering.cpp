/*
 * gto_ordering.cpp
 *
 *  Created on: Apr 27, 2020
 *      Author: JG
 */
#include "BasisSet.h"
#include "Shell.h"
#include "utils.h"
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include <tuple>
#include "gto_ordering.h"
#include <iostream>
#include <Eigen/Dense>


std::vector<std::array<size_t,3>> opencap_carts_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {{0,0,0}};
	//p
	if(shell.l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (shell.l==2)
		return {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
	//f
	else if (shell.l==3)
		return {{3,0,0},{2,1,0},{2,0,1},{1,0,2},{1,1,1},{1,2,0},{0,3,0},
				{0,2,1},{0,1,2},{0,0,3}};
	//g
	else if (shell.l==4)
		return { {4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},
		         {1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};
	else
	{
		return {{0,0,0}};
	}
}

std::vector<std::array<size_t,3>> pyscf_carts_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {{0,0,0}};
	//p
	if(shell.l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (shell.l==2)
		return {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
	//f
	else if (shell.l==3)
		return {{3,0,0},{2,1,0},{2,0,1},{1,0,2},{1,1,1},{1,2,0},{0,3,0},
				{0,2,1},{0,1,2},{0,0,3}};
	//g
	else if (shell.l==4)
		return { {4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},
		         {1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};
	else
	{
		return {{0,0,0}};
	}
}

std::vector<int> pyscf_harmonic_ordering(Shell shell)
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
	//h
	else if(shell.l==5)
		return {-5,-4,-3,-2,-1,0,1,2,3,4,5};
	else
		return {0,0,0};
}

std::vector<int> opencap_harmonic_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {0};
	//p
	if(shell.l==1)
		return {-1,0,1};
	//d
	else if (shell.l==2)
		return {-2,-1,0,1,2};
	//f
	else if (shell.l==3)
		return {-3,-2,-1,0,1,2,3};
	//g
	else if (shell.l==4)
		return {-4,-3,-2,-1,0,1,2,3,4};
	//h
	else if(shell.l==5)
		return {-5,-4,-3,-2,-1,0,1,2,3,4,5};
	else
		return {0,0,0};
}

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

// matrix in opencap ordering --> matrix in molden ordering
void to_molden_ordering(Eigen::MatrixXd &opencap_mat, BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> opencap_order = opencap_harmonic_ordering(shell);
			std::vector<int> molden_order = molden_harmonic_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(opencap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(opencap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	Eigen::MatrixXd per_mat(bs.Nbasis,bs.Nbasis);
	per_mat= Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	opencap_mat = per_mat.transpose()* opencap_mat * per_mat;
}

// matrix in opencap ordering --> matrix in qchem ordering
void to_qchem_ordering(Eigen::MatrixXd &opencap_mat, BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> opencap_order = opencap_harmonic_ordering(shell);
			std::vector<int> molden_order = qchem_harmonic_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(opencap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(opencap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	Eigen::MatrixXd per_mat(bs.Nbasis,bs.Nbasis);
	per_mat= Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	opencap_mat = per_mat.transpose()* opencap_mat * per_mat;
}

// matrix in opencap ordering --> matrix in qchem ordering
void molden_to_qchem_ordering(Eigen::MatrixXd &opencap_mat, BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> opencap_order = molden_harmonic_ordering(shell);
			std::vector<int> molden_order = qchem_harmonic_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(opencap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<molden_order.size();j++)
				{
					if(opencap_order[i]==molden_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	Eigen::MatrixXd per_mat(bs.Nbasis,bs.Nbasis);
	per_mat= Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	opencap_mat = per_mat.transpose()* opencap_mat * per_mat;
}

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

bool same_atom(std::array<double,3> shell_origin,std::array<double,3> atom_coords)
{
	for(size_t i=0;i<3;i++)
	{
		if(shell_origin[i]!=atom_coords[i])
			return false;
	}
	return true;
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
		while (!same_atom(shell.origin,geometry[atm_idx].coords))
			atm_idx++;
		reordered_shells[atm_idx][shell.l].push_back(shell);
	}
	return reordered_shells;
}

size_t find_matching_index(Shell shell, int angmom, BasisSet original_bs)
{
	std::vector<int> opencap_order = opencap_harmonic_ordering(shell);
	size_t matching_idx = 0;
	for (auto shell2:original_bs.basis)
	{
		if (shell==shell2)
			break;
		else
			matching_idx +=shell2.num_bf;
	}
	for (size_t i=0;i<opencap_order.size();i++)
	{
		if (opencap_order[i]==angmom)
			return matching_idx + i;
	}
	std::cout << "Something's gone wrong." << std::endl;
	return -1;
}


// matrix in opencap ordering --> matrix in molcas ordering
void to_molcas_ordering(Eigen::MatrixXd &opencap_mat, BasisSet bs, std::vector<Atom> geometry)
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
	Eigen::MatrixXd per_mat(bs.Nbasis,bs.Nbasis);
	per_mat= Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	opencap_mat = per_mat.transpose() * opencap_mat * per_mat;
}

void to_pyscf_ordering(Eigen::MatrixXd &opencap_mat,BasisSet bs)
{
	std::vector<std::tuple<int,int>> swap_indices;
	int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		if (shell.pure)
		{
			std::vector<int> opencap_order = opencap_harmonic_ordering(shell);
			std::vector<int> pyscf_order = pyscf_harmonic_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<pyscf_order.size();j++)
				{
					if(opencap_order[i]==pyscf_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(shell);
			std::vector<std::array<size_t,3>> pyscf_order = pyscf_carts_ordering(shell);
			for (size_t i=0;i<opencap_order.size();i++)
			{
				for(size_t j=0; j<pyscf_order.size();j++)
				{
					if(opencap_order[i]==pyscf_order[j])
						swap_indices.push_back(std::make_tuple(i+bf_idx, j+bf_idx));
				}
			}
		}
		bf_idx+=shell.num_bf;
	}
	Eigen::MatrixXd per_mat(bs.Nbasis,bs.Nbasis);
	per_mat= Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	opencap_mat = per_mat.transpose()* opencap_mat * per_mat;
}
