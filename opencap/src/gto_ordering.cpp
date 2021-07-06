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
 * gto_ordering.cpp
 */

#include "gto_ordering.h"

#include <h5pp/h5pp.h>
#include <math.h>
#include <algorithm>
#include <array>
#include <Eigen/Dense>
#include <iostream>
#include <tuple>
#include <vector>

#include "BasisSet.h"
#include "opencap_exception.h"
#include "Shell.h"
#include "utils.h"


//from:http://cheminf.cmbi.ru.nl/molden/molden_format.html
std::vector<int> opencap_harmonic_ordering(int l)
{
	//s
	if(l == 0)
		return {0};
	//p
	if(l==1)
		return {1,-1,0};
	//d
	else if (l==2)
		return {0,1,-1,2,-2};
	//f
	else if (l==3)
		return {0,1,-1,2,-2,3,-3};
	//g
	else if (l==4)
		return {0,1,-1,2,-2,3,-3,4,-4};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

//from: http://cheminf.cmbi.ru.nl/molden/molden_format.html
std::vector<std::array<size_t,3>> opencap_carts_ordering(int l)
{
	//s
	if(l == 0)
		return {{0,0,0}};
	//p
	if(l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (l==2)
		return {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
	//f
	else if (l==3)
		return {{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},{2,0,1},{1,0,2},
				{0,1,2},{0,2,1},{1,1,1}};
	//g
	else if(l==4)
		return {{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},{1,3,0},{0,3,1},
		        {1,0,3},{0,1,3},{2,2,0},{2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

std::vector<int> psi4_harmonic_ordering(int l)
{
	//s
	if(l == 0)
		return {0};
	//p
	if(l==1)
		return {0,1,-1};
	//d
	else if (l==2)
		return {0,1,-1,2,-2};
	//f
	else if (l==3)
		return {0,1,-1,2,-2,3,-3};
	//g
	else if (l==4)
		return {0,1,-1,2,-2,3,-3,4,-4};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

// From: psi4/psi4/src/psi4/libmints/cartesianiter.cc
std::vector<std::array<size_t,3>> psi4_carts_ordering(int l)
{
	//s
	if(l == 0)
		return {{0,0,0}};
	//p
	if(l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (l==2)
		return {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
	//f
	else if (l==3)
		return {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},
			    {1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}};
	//g
	else if(l==4)
		return {{4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},
			    {2,0,2},{1,3,0},{1,2,1},{1,1,2},{1,0,3},
				{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}


// From: pyscf tools/molden.py
std::vector<std::array<size_t,3>> pyscf_carts_ordering(int l)
{
	//s
	if(l == 0)
		return {{0,0,0}};
	//p
	if(l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (l==2)
		return {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
	//f
	else if (l==3)
		return {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},{0,3,0},
				{0,2,1},{0,1,2},{0,0,3}};
	//g
	else if(l==4)
		return {{4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},
		        {1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

std::vector<int> pyscf_harmonic_ordering(int l)
{
	//s
	if(l == 0)
		return {0};
	//p
	if(l==1)
		return {1,-1,0};
	//d
	else if (l==2)
		return {-2,-1,0,1,2};
	//f
	else if (l==3)
		return {-3,-2,-1,0,1,2,3};
	//g
	else if (l==4)
		return {-4,-3,-2,-1,0,1,2,3,4};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

std::vector<bf_id> get_pyscf_ids(BasisSet &bs)
{
	std::vector<bf_id> ids;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		if (my_shell.pure)
		{
			std::vector<int> pyscf_order = pyscf_harmonic_ordering(my_shell.l);
			for(int m:pyscf_order)
				ids.push_back(bf_id(bs.shell_ids[i],m));
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(my_shell.l);
			std::vector<std::array<size_t,3>> pyscf_order = pyscf_carts_ordering(my_shell.l);
			for(auto m:pyscf_order)
			{
				auto it = std::find(opencap_order.begin(), opencap_order.end(), m);
				if (it != opencap_order.end())
				{
					int index = std::distance(opencap_order.begin(), it);
					ids.push_back(bf_id(bs.shell_ids[i],index));
				}
				else
				    opencap_throw("Something's gone wrong.");
			}
		}
	}
	return ids;
}

//same as molden except for g-type orbitals
std::vector<std::array<size_t,3>> qchem_carts_ordering(int l)
{
	//s
	if(l == 0)
		return {{0,0,0}};
	//p
	if(l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (l==2)
		return {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
	//f
	else if (l==3)
		return {{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},{2,0,1},{1,0,2},
				{0,1,2},{0,2,1},{1,1,1}};
	//g
	else if(l==4)
		return {{0,0,4},{0,1,3},{0,2,2},{0,3,1},{0,4,0},{1,0,3},{1,1,2},
		        {1,2,1},{1,3,0},{2,0,2},{2,1,1},{2,2,0},{3,0,1},{3,1,0},{4,0,0}};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

//should be same as molden
std::vector<int> qchem_harmonic_ordering(int l)
{
	//s
	if(l == 0)
		return {0};
	//p
	if(l==1)
		return {1,-1,0};
	//d
	else if (l==2)
		return {0,1,-1,2,-2};
	//f
	else if (l==3)
		return {0,1,-1,2,-2,3,-3};
	//g
	else if (l==4)
		return {0,1,-1,2,-2,3,-3,4,-4};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

std::vector<bf_id> get_qchem_ids(BasisSet &bs)
{
	std::vector<bf_id> ids;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		if (my_shell.pure)
		{
			std::vector<int> qchem_order = qchem_harmonic_ordering(my_shell.l);
			for(int m:qchem_order)
				ids.push_back(bf_id(bs.shell_ids[i],m));
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(my_shell.l);
			std::vector<std::array<size_t,3>> qchem_order = qchem_carts_ordering(my_shell.l);
			for(auto m:qchem_order)
			{
				auto it = std::find(opencap_order.begin(), opencap_order.end(), m);
				if (it != opencap_order.end())
				{
					int index = std::distance(opencap_order.begin(), it);
					ids.push_back(bf_id(bs.shell_ids[i],index));
				}
				else
				    opencap_throw("Something's gone wrong.");
			}
		}
	}
	return ids;
}

std::vector<bf_id> get_psi4_ids(BasisSet &bs)
{
	std::vector<bf_id> ids;
	for(size_t i=0;i<bs.basis.size();i++)
	{
		Shell my_shell = bs.basis[i];
		if (my_shell.pure)
		{
			std::vector<int> qchem_order = psi4_harmonic_ordering(my_shell.l);
			for(int m:qchem_order)
				ids.push_back(bf_id(bs.shell_ids[i],m));
		}
		else
		{
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(my_shell.l);
			std::vector<std::array<size_t,3>> qchem_order = psi4_carts_ordering(my_shell.l);
			for(auto m:qchem_order)
			{
				auto it = std::find(opencap_order.begin(), opencap_order.end(), m);
				if (it != opencap_order.end())
				{
					int index = std::distance(opencap_order.begin(), it);
					ids.push_back(bf_id(bs.shell_ids[i],index));
				}
				else
				    opencap_throw("Something's gone wrong.");
			}
		}
	}
	return ids;
}

std::vector<std::array<size_t,3>> molcas_carts_ordering(int l)
{
	//s
	if(l == 0)
		return {{0,0,0}};
	//p
	if(l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (l==2)
		return {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
	//f
	else if (l==3)
		return {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},{0,3,0},
				{0,2,1},{0,1,2},{0,0,3}};
	//g
	else if(l==4)
		return { {4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},
		         {1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

std::vector<int> molcas_harmonic_ordering(int l)
{
	//s
	if(l == 0)
		return {0};
	//p
	if(l==1)
		return {1,-1,0};
	//d
	else if (l==2)
		return {-2,-1,0,1,2};
	//f
	else if (l==3)
		return {-3,-2,-1,0,1,2,3};
	//g
	else if (l==4)
		return {-4,-3,-2,-1,0,1,2,3,4};
	else
		opencap_throw("Error: Only up to G type orbitals are supported.");
}

size_t get_molcas_cart_idx(int l,int m)
{
	int cur = -1*(opencap_harmonic_ordering(l).size()/2);
	size_t num_functions = molcas_carts_ordering(l).size();
	std::vector<int> indices;
	for(size_t i=1;i<=num_functions;i++)
	{
		indices.push_back(cur);
		cur++;
	}
	auto it = std::find(indices.begin(), indices.end(), m);
	if (it != indices.end())
	{
		return std::distance(indices.begin(), it);
	}
	else
	    opencap_throw("Something's gone wrong.");
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
		int l = basis_ids(i,2);
		size_t shell_idx = bs.get_index_of_shell_id(shell_id(ctr,shell_num,l));
		Shell my_shell = bs.basis[shell_idx];
		long m = basis_ids(i,3);
		//cartesian, need to translate molcas cartesian IDs to opencap IDs
		if(l<0)
		{
			std::vector<std::array<size_t,3>> molcas_order = molcas_carts_ordering(my_shell.l);
			std::vector<std::array<size_t,3>> opencap_order = opencap_carts_ordering(my_shell.l);
			std::array<size_t,3> cart_m = molcas_order[get_molcas_cart_idx(my_shell.l,m)];
			auto it = std::find(opencap_order.begin(), opencap_order.end(), cart_m);
			if (it != opencap_order.end())
			{
				int index = std::distance(opencap_order.begin(), it);
				ids.push_back(bf_id(bs.shell_ids[shell_idx],index));
			}
			else
			    opencap_throw("I didn't find that shell.");
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

void reorder_matrix(Eigen::MatrixXd &mat,std::vector<bf_id> &original_ids,std::vector<bf_id> &final_ids)
{
    std::vector<std::tuple<int,int>> swap_indices;
    for(size_t i=0;i<original_ids.size();i++)
    {
        for(size_t j=0;j<original_ids.size();j++)
        {
            if(original_ids[i]==final_ids[j])
            {
                swap_indices.push_back(std::make_tuple(i,j));
            }
        }
    }
    Eigen::MatrixXd per_mat(original_ids.size(),original_ids.size());
    per_mat= Eigen::MatrixXd::Zero(original_ids.size(),original_ids.size());
    for(auto t:swap_indices)
    per_mat(std::get<0>(t),std::get<1>(t))=1;
    // permute indices: P^T * A * P
    mat = per_mat.transpose()* mat * per_mat;
}
