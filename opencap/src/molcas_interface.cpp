/*
 * readMolcasHDF5.cpp
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */
#include <algorithm>
#include <iostream>
#include <fstream>
#include "molcas_interface.h"
#include "opencap_exception.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <h5pp/h5pp.h>
#include "gto_ordering.h"
#include "BasisSet.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar,int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}

std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> read_rassi_tdms(std::string filename,BasisSet bs)
{
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	//first lets check if the dimensions are correct
	std::vector<long> nbas_vec = file.readAttribute<std::vector<long>>("NBAS","/");
	long nbas=0;
	for(auto num_bas:nbas_vec)
		nbas+=num_bas;
	if(nbas!=bs.Nbasis)
		opencap_throw("Error: dimensions of TDMs do not match specified basis set.");
	//now lets load the densities
	Eigen::Tensor<double,3> rass_data, spin_dens;
	std::vector<std::vector<Eigen::MatrixXd>> alpha_opdms;
	std::vector<std::vector<Eigen::MatrixXd>> beta_opdms;
	file.readDataset(rass_data,"SFS_TRANSITION_DENSITIES");
	file.readDataset(spin_dens,"SFS_TRANSITION_SPIN_DENSITIES");
	//now lets loop over the matrices
    const auto& d = rass_data.dimensions();
    for(long i=0;i<d[0];i++)
    {
		std::vector<Eigen::MatrixXd> alpha_state_row;
		std::vector<Eigen::MatrixXd> beta_state_row;
		for (long j=0;j<d[1];j++)
		{
			long nbas = sqrt(d[2]);
			Eigen::array<long,3> offset = {i,j,0};    //Starting point
			Eigen::array<long,3> extent = {1,1,nbas}; //end point
			Eigen::Tensor<double, 2> dmt_slice = rass_data.slice(offset, extent).reshape(Eigen::array<long,2>{nbas,nbas});
			Eigen::MatrixXd  dmt_mat =  Tensor_to_Matrix(dmt_slice,nbas,nbas);
			Eigen::Tensor<double, 2> spin_slice = spin_dens.slice(offset, extent).reshape(Eigen::array<long,2>{nbas,nbas});
			Eigen::MatrixXd  spin_mat =  Tensor_to_Matrix(dmt_slice,nbas,nbas);
			Eigen::MatrixXd alpha_opdm,beta_opdm;
			alpha_opdm= 0.5*(dmt_mat+spin_mat);
			beta_opdm = 0.5*(dmt_mat-spin_mat);
			to_opencap_ordering(alpha_opdm,bs,get_molcas_ids(bs,filename));
			to_opencap_ordering(beta_opdm,bs,get_molcas_ids(bs,filename));
			alpha_state_row.push_back(alpha_opdm);
			beta_state_row.push_back(beta_opdm);
		}
		alpha_opdms.push_back(alpha_state_row);
		beta_opdms.push_back(beta_state_row);
    }
    //symmetric, for some reason molcas does this wrong
    for (size_t i=0;i<d[0];i++)
    {
    	for(size_t j=0;j<i;j++)
    	{
    		alpha_opdms[i][j]= alpha_opdms[j][i];
    		beta_opdms[i][j]= beta_opdms[j][i];
    	}
    }
    return {alpha_opdms,beta_opdms};
}

Eigen::Map<const Eigen::MatrixXd> reshape (const Eigen::VectorXd& b, const uint n, const uint m) {
    return Eigen::Map<const Eigen::MatrixXd>(b.data(), n, m);
}

Eigen::MatrixXd read_rassi_overlap(std::string filename,BasisSet bs)
{
	Eigen::VectorXd overlap_vec;
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	file.readDataset(overlap_vec,"AO_OVERLAP_MATRIX");
	auto nsym = file.readAttribute<long>("NSYM", "/");
	if(nsym >1)
	{
		Eigen::VectorXd desym_vec;
		file.readDataset(desym_vec,"DESYM_MATRIX");
		Eigen::MatrixXd desym_mat;
		std::vector<long> nbas_vec = file.readAttribute<std::vector<long>>("NBAS","/");
		long nbas=0;
		for(auto num_bas:nbas_vec)
			nbas+=num_bas;
		desym_mat = reshape(desym_vec,nbas,nbas);
		Eigen::MatrixXd smat;
		smat = Eigen::MatrixXd::Zero(nbas,nbas);
		//fill the blocks
		size_t elements_index = 0;
		size_t bf_index = 0;
		for(size_t isym=0;isym<nsym;isym++)
		{
			size_t n_elements = nbas_vec[isym] * nbas_vec[isym];
			Eigen::VectorXd block_vec(n_elements);
			Eigen::MatrixXd block_mat;
			//grab the elements we need
			for(size_t i=0;i<n_elements;i++)
				block_vec(i) = overlap_vec(elements_index+i);
			block_mat = reshape(block_vec,nbas_vec[isym],nbas_vec[isym]);
			smat.block(bf_index,bf_index,nbas_vec[isym],nbas_vec[isym]) = block_mat;
			//update our index counters
			elements_index+=n_elements;
			bf_index+= nbas_vec[isym];
		}
		Eigen::MatrixXd desym_overlap = desym_mat * smat * desym_mat.transpose();
		to_opencap_ordering(desym_overlap,bs,get_molcas_ids(bs,filename));
		return desym_overlap;
	}
	else
	{
		int nbas = sqrt(overlap_vec.size());
		Eigen::MatrixXd smat;
		smat = reshape(overlap_vec,nbas,nbas);
		to_opencap_ordering(smat,bs,get_molcas_ids(bs,filename));
		return smat;
	}
}

Eigen::MatrixXd read_mscaspt2_heff(size_t nstates, std::string filename)
{
	Eigen::MatrixXd ZERO_ORDER_H(nstates,nstates);
	ZERO_ORDER_H= Eigen::MatrixXd::Zero(nstates,nstates);
	std::ifstream is(filename);
	if (is.good())
	{
		std::string line, rest;
		std::getline(is,line);
		while (line.find("MULTI-STATE CASPT2 SECTION")== std::string::npos)
			std::getline(is,line);
		//get diagonal shift
		for (size_t i=1;i<=3;i++)
			std::getline(is,line);
		std::vector<std::string> split_line = split(line,' ');
		double E_shift = std::stod(split_line[split_line.size()-1]);
		for(size_t i=1;i<=2;i++)
			std::getline(is,line);
		for (size_t i=1;i<=(nstates + 5 - 1) / 5;i++)
		{
			for (size_t j=1;j<=2;j++)
				std::getline(is,line);
			//now time to start reading in the matrix elements
			for (size_t j=1;j<=nstates-((i-1)*5);j++)
			{
				std::getline(is,line);
				std::vector<std::string> tokens = split(line,' ');
				size_t row_idx = std::stoul(tokens[0]);
				for(size_t k=1;k<tokens.size();k++)
				{
					ZERO_ORDER_H(row_idx-1,k-1+(i-1)*5)=std::stod(tokens[k]);
					ZERO_ORDER_H(k-1+(i-1)*5,row_idx-1)=std::stod(tokens[k]);
				}
			}
		}
		for (size_t i=0;i<nstates;i++)
			ZERO_ORDER_H(i,i)+=E_shift;
	}
	return ZERO_ORDER_H;
}

std::vector<Atom> read_geometry_from_rassi(std::string filename)
{
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	std::vector<Atom> atoms;
	std::vector<long> atomic_nums;
	Eigen::MatrixXd coords_mat;
	file.readDataset(atomic_nums,"CENTER_ATNUMS");
	file.readDataset(coords_mat,"CENTER_COORDINATES");
	for(size_t i=0;i<coords_mat.rows();i++)
		atoms.push_back(Atom(atomic_nums[i],coords_mat(i,0),coords_mat(i,1),coords_mat(i,2)));
	return atoms;
}

BasisSet read_basis_from_rassi(std::string filename,std::vector<Atom> atoms)
{
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	Eigen::Matrix<long,Eigen::Dynamic,Eigen::Dynamic> prim_ids,basis_ids;
	Eigen::MatrixXd prims;
	auto nsym = file.readAttribute<long>("NSYM", "/");
	file.readDataset(prims,"PRIMITIVES");
	file.readDataset(prim_ids,"PRIMITIVE_IDS");
	std::string bf_ids_tag;
	if(nsym>1)
		bf_ids_tag="DESYM_BASIS_FUNCTION_IDS";
	else
		bf_ids_tag ="BASIS_FUNCTION_IDS";
	file.readDataset(basis_ids,bf_ids_tag);
	BasisSet bs;
	for(auto atm:atoms)
		bs.centers.push_back(atm.coords);
	for(size_t i=0;i<basis_ids.rows();i++)
	{
		long ctr = basis_ids(i,0);
		long shell_num = basis_ids(i,1);
		long l = basis_ids(i,2);
		long m = basis_ids(i,3);
		shell_id id(ctr,shell_num,l);
		int bs_idx = bs.get_index_of_shell_id(id);
		if (bs_idx==-1)
		{
			Shell new_shell(l,bs.centers[ctr-1]);
			new_shell.l=abs(l);
			if(l<0 && abs(l)>1)
				new_shell.pure=false;
			bs.add_shell(new_shell);
		}
	}
	for(size_t i=0;i<prim_ids.rows();i++)
	{
		long ctr = prim_ids(i,0);
		long l = prim_ids(i,1);
		long shell_num = prim_ids(i,2);
		double exp = prims(i,0);
		double coeff = prims(i,1);
		shell_id id(ctr,shell_num,l);
		int bs_idx = bs.get_index_of_shell_id(id);
		if (coeff!=0)
			bs.basis[bs_idx].add_primitive(exp,coeff);
	}
	bs.normalize();
	return bs;
}
