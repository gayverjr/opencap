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
#include "opencap_exception.h"
#include <iostream>
#include <Eigen/Dense>


std::vector<std::array<size_t,3>> opencap_carts_ordering(Shell &shell)
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

std::vector<int> opencap_harmonic_ordering(Shell &shell)
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

std::vector<std::array<size_t,3>> pyscf_carts_ordering(Shell &shell)
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

std::vector<int> pyscf_harmonic_ordering(Shell &shell)
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

std::vector<bf_id> get_pyscf_ids(BasisSet &bs)
{
	std::vector<bf_id> ids;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		if (my_shell.pure)
		{
			std::vector<int> pyscf_order = pyscf_harmonic_ordering(my_shell);
			for(int m:pyscf_order)
				ids.push_back(bf_id(bs.shell_ids[i],m));
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(my_shell);
			std::vector<std::array<size_t,3>> pyscf_order = pyscf_carts_ordering(my_shell);
			for(auto m:pyscf_order)
			{
				auto it = std::find(opencap_order.begin(), opencap_order.end(), m);
				if (it != opencap_order.end())
				{
					int index = std::distance(opencap_order.begin(), it);
					int m_id = index-(opencap_order.size()-1)/2;
					ids.push_back(bf_id(bs.shell_ids[i],m_id));
				}
				else
				    opencap_throw("Something's gone wrong.");
			}
		}
	}
	return ids;
}

std::vector<std::array<size_t,3>> molden_carts_ordering(Shell &shell)
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

std::vector<int> molden_harmonic_ordering(Shell &shell)
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

std::vector<bf_id> get_molden_ids(BasisSet &bs)
{
	std::vector<bf_id> ids;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		if (my_shell.pure)
		{
			std::vector<int> molden_order = molden_harmonic_ordering(my_shell);
			for(int m:molden_order)
				ids.push_back(bf_id(bs.shell_ids[i],m));
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(my_shell);
			std::vector<std::array<size_t,3>> molden_order = molden_carts_ordering(my_shell);
			for(auto m:molden_order)
			{
				auto it = std::find(opencap_order.begin(), opencap_order.end(), m);
				if (it != opencap_order.end())
				{
					int index = std::distance(opencap_order.begin(), it);
					int m_id = index-(opencap_order.size()-1)/2;
					ids.push_back(bf_id(bs.shell_ids[i],m_id));
				}
				else
				    opencap_throw("Something's gone wrong.");
			}
		}
	}
	return ids;
}

std::vector<std::array<size_t,3>> molcas_carts_ordering(Shell &shell)
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

std::vector<int> molcas_harmonic_ordering(Shell &shell)
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

std::vector<bf_id> get_molcas_ids(BasisSet &bs,std::string rassi_filename)
{
	std::vector<bf_id> ids;
	h5pp::File file(rassi_filename, h5pp::FilePermission::READONLY);
	Eigen::Matrix<long,Eigen::Dynamic,Eigen::Dynamic> basis_ids;
	auto nsym = file.readAttribute<long>("NSYM", "/");
	std::string bf_ids_tag;
	if(nsym>1)
		bf_ids_tag="DESYM_BASIS_FUNCTION_IDS";
	else
		bf_ids_tag ="BASIS_FUNCTION_IDS";
	file.readDataset(basis_ids,bf_ids_tag);
	for(size_t i=0;i<basis_ids.rows();i++)
	{
		long ctr = basis_ids(i,0);
		long shell_num = basis_ids(i,1);
		long l = basis_ids(i,2);
		size_t shell_idx = bs.get_index_of_shell_id(shell_id(ctr,shell_num,l));
		Shell my_shell = bs.basis[shell_idx];
		long m = basis_ids(i,3);
		//cartesian, need to translate molcas cartesian IDs to opencap IDs
		if(l<0)
		{
			std::vector<std::array<size_t,3>> molcas_order = molcas_carts_ordering(my_shell);
			std::vector<std::array<size_t,3>> opencap_order = molcas_carts_ordering(my_shell);
			std::array<size_t,3> cart_m = molcas_order[m+(molcas_order.size()-1)/2];
			auto it = std::find(opencap_order.begin(), opencap_order.end(), cart_m);
			if (it != opencap_order.end())
			{
				int index = std::distance(opencap_order.begin(), it);
				int m_id = index-(opencap_order.size()-1)/2;
				ids.push_back(bf_id(bs.shell_ids[i],m_id));
			}
			else
			    opencap_throw("Something's gone wrong.");
		}
		//spherical harmonic, don't need to do anything fancy
		else
		{
			ids.push_back(bf_id(bs.shell_ids[shell_idx],m));
		}
	}
	return ids;
}

void to_opencap_ordering(Eigen::MatrixXd &mat,BasisSet &bs,std::vector<bf_id> input_ids)
{
	std::vector<std::tuple<int,int>> swap_indices;
	for(size_t i=0;i<input_ids.size();i++)
	{
		for(size_t j=0;j<bs.bf_ids.size();j++)
		{
			if(input_ids[i]==bs.bf_ids[j])
			{
				swap_indices.push_back(std::make_tuple(i, j));
			}
		}
	}
	Eigen::MatrixXd per_mat(bs.Nbasis,bs.Nbasis);
	per_mat= Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
	for(auto t:swap_indices)
		per_mat(std::get<0>(t),std::get<1>(t))=1;
	// permute indices: P^T * A * P
	mat = per_mat.transpose()* mat * per_mat;
}
