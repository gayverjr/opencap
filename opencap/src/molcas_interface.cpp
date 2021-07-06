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
 * molcas_interface.cpp
 */

#include "molcas_interface.h"

#include <h5pp/h5pp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>

#include "BasisSet.h"
#include "gto_ordering.h"
#include "opencap_exception.h"
#include "utils.h"

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar,int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}

Eigen::Map<const Eigen::MatrixXd> reshape (const Eigen::VectorXd& b, const uint n, const uint m) {
    return Eigen::Map<const Eigen::MatrixXd>(b.data(), n, m);
}

void read_rassi_tdms(std::vector<std::vector<Eigen::MatrixXd>> &alpha_opdms,
		std::vector<std::vector<Eigen::MatrixXd>> &beta_opdms,
		std::string filename, BasisSet bs,size_t nstates)
{
    alpha_opdms = std::vector< std::vector<Eigen::MatrixXd>>(nstates, std::vector<Eigen::MatrixXd> (nstates));
    beta_opdms = std::vector< std::vector<Eigen::MatrixXd>>(nstates, std::vector<Eigen::MatrixXd> (nstates));
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	//first lets check if the dimensions are correct
	std::vector<long> nbas_vec = file.readAttribute<std::vector<long>>("NBAS","/");
	long nbas=0;
	for(auto num_bas:nbas_vec)
		nbas+=num_bas;
	if(nbas!=bs.Nbasis)
		opencap_throw("Error: dimensions of RASSI basis set does not match specified basis set.");
	auto nsym = file.readAttribute<long>("NSYM", "/");
	//now lets load the densities
	Eigen::Tensor<double,3> rass_data, spin_dens;
	try
	{
		file.readDataset(rass_data,"SFS_TRANSITION_DENSITIES");
		file.readDataset(spin_dens,"SFS_TRANSITION_SPIN_DENSITIES");
	}
	catch(exception &e)
	{
		opencap_throw("Error: SFS_TRANSITION_DENSITIES dataset not found. Use the TRD1 keyword"
				" in the RASSI module to activate transition densities.")
	}
    const auto& d = rass_data.dimensions();

    if (nsym>1)
    {
    	//Step 1: get the desymmetrization matrix
		Eigen::VectorXd desym_vec;
		file.readDataset(desym_vec,"DESYM_MATRIX");
		Eigen::MatrixXd desym_mat;
		std::vector<long> nbas_vec = file.readAttribute<std::vector<long>>("NBAS","/");
		desym_mat = reshape(desym_vec,nbas,nbas);
		long mat_size=0;
		for(auto num_bas:nbas_vec)
			mat_size+=num_bas*num_bas;
		//now lets loop over the matrices
		for(long i=0;i<d[0];i++)
		{
            if(d[0]!=nstates)
                opencap_throw("Error: Found " + std::to_string(d[0]) + " states in RASSI file, but "
                          + std::to_string(nstates) + " states were specified.");
			std::vector<Eigen::MatrixXd> alpha_state_row;
			std::vector<Eigen::MatrixXd> beta_state_row;
			for (long j=0;j<d[1];j++)
			{
				//step 2: get the raw matrices
				Eigen::array<long,3> offset = {i,j,0};    //Starting point
				Eigen::array<long,3> extent = {1,1,1}; //end point
				Eigen::Tensor<double, 2> dmt_slice = rass_data.slice(offset, extent).reshape(Eigen::array<long,2>{mat_size,(long)1});
				Eigen::MatrixXd  dmt_vec =  Tensor_to_Matrix(dmt_slice,mat_size,(long)1);
				Eigen::Tensor<double, 2> spin_slice = spin_dens.slice(offset, extent).reshape(Eigen::array<long,2>{mat_size,(long)1});
				Eigen::MatrixXd  spin_vec =  Tensor_to_Matrix(dmt_slice,mat_size,(long)1);

				//step 3: turn into bigger matrix by adding zeros
				Eigen::MatrixXd dmat,spin;
				dmat = Eigen::MatrixXd::Zero(nbas,nbas);
				spin = Eigen::MatrixXd::Zero(nbas,nbas);
				//fill the blocks
				size_t elements_index = 0;
				size_t bf_index = 0;
				for(size_t isym=0;isym<nsym;isym++)
				{
					size_t n_elements = nbas_vec[isym] * nbas_vec[isym];
					Eigen::VectorXd dmat_block_vec(n_elements),spin_block_vec(n_elements);
					Eigen::MatrixXd dmat_block,spin_block;
					//grab the elements we need
					for(size_t k=0;k<n_elements;k++)
					{
						dmat_block_vec(k) = dmt_vec(elements_index+k);
						spin_block_vec(k) = spin_vec(elements_index+k);
					}
					dmat_block = reshape(dmat_block_vec,nbas_vec[isym],nbas_vec[isym]);
					spin_block = reshape(spin_block_vec,nbas_vec[isym],nbas_vec[isym]);
					dmat.block(bf_index,bf_index,nbas_vec[isym],nbas_vec[isym]) = dmat_block;
					spin.block(bf_index,bf_index,nbas_vec[isym],nbas_vec[isym]) = spin_block;
					//update our index counters
					elements_index+=n_elements;
					bf_index+= nbas_vec[isym];
				}
                
				//step 4: desymmetrize
				Eigen::MatrixXd desym_dmat = desym_mat * dmat * desym_mat.transpose();
				Eigen::MatrixXd desym_spin = desym_mat * spin * desym_mat.transpose();

				//step 5: save
				Eigen::MatrixXd alpha_opdm,beta_opdm;
				alpha_opdm= 0.5*(desym_dmat+desym_spin);
				beta_opdm = 0.5*(desym_dmat-desym_spin);
				to_opencap_ordering(alpha_opdm,bs,get_molcas_ids(bs,filename));
				to_opencap_ordering(beta_opdm,bs,get_molcas_ids(bs,filename));
                alpha_opdms[i][j] = alpha_opdm;
                beta_opdms[i][j] = beta_opdm;
			}
		}

    }
    else
    {
		//now lets loop over the matrices
		for(long i=0;i<d[0];i++)
		{
            if(d[0]!=nstates)
            opencap_throw("Error: Found " + std::to_string(d[0]) + " states in RASSI file, but "
                          + std::to_string(nstates) + " states were specified.");
			for (long j=0;j<d[1];j++)
			{
				Eigen::array<long,3> offset = {i,j,0};    //Starting point
				Eigen::array<long,3> extent = {1,1,1}; //end point
				Eigen::Tensor<double, 2> dmt_slice = rass_data.slice(offset, extent).reshape(Eigen::array<long,2>{nbas,nbas});
				Eigen::MatrixXd  dmt_mat =  Tensor_to_Matrix(dmt_slice,nbas,nbas);
				Eigen::Tensor<double, 2> spin_slice = spin_dens.slice(offset, extent).reshape(Eigen::array<long,2>{nbas,nbas});
				Eigen::MatrixXd  spin_mat =  Tensor_to_Matrix(dmt_slice,nbas,nbas);
				Eigen::MatrixXd alpha_opdm,beta_opdm;
				alpha_opdm= 0.5*(dmt_mat+spin_mat);
				beta_opdm = 0.5*(dmt_mat-spin_mat);
				to_opencap_ordering(alpha_opdm,bs,get_molcas_ids(bs,filename));
				to_opencap_ordering(beta_opdm,bs,get_molcas_ids(bs,filename));
                alpha_opdms[i][j] = alpha_opdm;
                beta_opdms[i][j] = beta_opdm;
			}
		}
    }

    std::cout << "Warning: TDM M-->N is assumed to be conjugate transpose of "
    << "TDM N-->M where M>N" << std::endl;
    for (size_t i=0;i<nstates;i++)
    {
    	for(size_t j=i+1;j<nstates;j++)
    	{
    		alpha_opdms[j][i]= alpha_opdms[i][j].adjoint();
    		beta_opdms[j][i]= beta_opdms[i][j].adjoint();
    	}
    }
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

Eigen::MatrixXd read_rotation_matrix(size_t nstates, std::ifstream &is)
{
	Eigen::MatrixXd rotation_matrix(nstates,nstates);
	rotation_matrix = Eigen::MatrixXd::Zero(nstates,nstates);
	std::string line, rest;
	size_t num_groups = nstates%5==0 ? nstates/5 : nstates/5+1;
	for (size_t i=1;i<=num_groups;i++)
	{
		std::getline(is,line);
        if(i>1)
            std::getline(is,line);
		for (size_t j=1;j<=nstates;j++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			size_t row_idx = std::stoul(tokens[0]);
			for(size_t k=1;k<tokens.size();k++)
			{
				size_t col_idx = k-1+(i-1)*5;
				if(col_idx>=rotation_matrix.cols() || row_idx-1 >= rotation_matrix.rows())
					opencap_throw("Error: State index of out bounds. There is a problem with the OpenMolcas output file. Exiting...");
				rotation_matrix(row_idx-1,col_idx)=std::stod(tokens[k]);
			}
		}
	}
	return rotation_matrix;
}

Eigen::MatrixXd read_mscaspt2_heff(size_t nstates, std::string filename)
{
	Eigen::MatrixXd ZERO_ORDER_H(nstates,nstates);
	Eigen::MatrixXd rotation_matrix(nstates,nstates);
	ZERO_ORDER_H= Eigen::MatrixXd::Zero(nstates,nstates);
	rotation_matrix = Eigen::MatrixXd::Zero(nstates,nstates);
	for (size_t i=0;i<nstates;i++)
		rotation_matrix(i,i)=1.0;
	std::ifstream is(filename);
	if (is.good())
	{
		std::string line, rest;
		std::getline(is,line);
		while (line.find("Number of CI roots used")== std::string::npos && is.peek()!=EOF)
			std::getline(is,line);
		size_t num_states = stoi(split(line,' ').back());
		if (num_states!=nstates)
			opencap_throw("Error: "+std::to_string(num_states)+ " roots were found in the OpenMolcas "
					"output file, but " + std::to_string(nstates) +" states were specified in the input. "
							"Exiting...");
		while (line.find("MULTI-STATE CASPT2 SECTION")== std::string::npos && is.peek()!=EOF)
		{
			std::getline(is,line);
			if (line.find("H0 eigenvectors:")!= std::string::npos)
			{
				std::cout << "Warning: rotating effective Hamiltonian into basis "
						  << "of original CASCI states." << std::endl;
				rotation_matrix = read_rotation_matrix(nstates,is);
			}
		}
		if (is.peek()==EOF)
			opencap_throw("Error: Reached end of file before MULTI-STATE CASPT2 SECTION.");
		//get diagonal shift
		for (size_t i=1;i<=3;i++)
			std::getline(is,line);
        double E_shift = 0.0;
        if(line.find("Output diagonal")!=std::string::npos)
        {
            std::vector<std::string> split_line = split(line,' ');
            E_shift = std::stod(split_line[split_line.size()-1]);
            for(size_t i=1;i<=2;i++)
                std::getline(is,line);
        }
		size_t num_groups = nstates%5==0 ? nstates/5 : nstates/5+1;
		for (size_t i=1;i<=num_groups;i++)
		{
			for (size_t j=1;j<=2;j++)
				std::getline(is,line);
			//now time to start reading in the matrix elements
			size_t states_in_group = nstates - (i-1)*5;
			for (size_t j=1;j<=states_in_group;j++)
			{
				std::getline(is,line);
				std::vector<std::string> tokens = split(line,' ');
				size_t row_idx = std::stoul(tokens[0]);
				for(size_t k=1;k<tokens.size();k++)
				{
					size_t col_idx = k-1+(i-1)*5;
					if(col_idx>=ZERO_ORDER_H.cols() || row_idx-1 >= ZERO_ORDER_H.rows())
						opencap_throw("Error: State index of out bounds. There is a problem with the OpenMolcas output file. Exiting...");
					ZERO_ORDER_H(row_idx-1,col_idx)=std::stod(tokens[k]);
					ZERO_ORDER_H(col_idx,row_idx-1)=std::stod(tokens[k]);
				}
			}
		}
		for (size_t i=0;i<nstates;i++)
			ZERO_ORDER_H(i,i)+=E_shift;
	}
	else
	{
    	opencap_throw("Error: I couldn't read:" + filename);
	}
	return rotation_matrix*ZERO_ORDER_H*rotation_matrix.transpose();
}

Eigen::MatrixXd read_nevpt2_heff(size_t nstates, std::string filename, std::string method)
{
	std::string line_to_find;
	if(compare_strings(method,"pc-nevpt2"))
		line_to_find = "Zero + second order effective Hamiltonian (PC)";
	else if(compare_strings(method,"sc-nevpt2"))
		line_to_find = "Zero + second order Effective Hamiltonian (SC)";
	else
		opencap_throw("Either pc-nevpt2 or sc-nevpt2 should be selected.");
	Eigen::MatrixXd ZERO_ORDER_H(nstates,nstates);
	ZERO_ORDER_H= Eigen::MatrixXd::Zero(nstates,nstates);
	std::ifstream is(filename);
	if (is.good())
	{
		std::string line, rest;
		std::getline(is,line);
		while (line.find(line_to_find)== std::string::npos && is.peek()!=EOF)
			std::getline(is,line);
		if(is.peek()==EOF)
			opencap_throw("Error: Unable to find QD-NEVPT2 effective Hamiltonian.");
		std::getline(is,line);
		std::getline(is,line);
		size_t num_groups = nstates%5==0 ? nstates/5 : nstates/5+1;
		for (size_t i=1;i<=num_groups;i++)
		{
			std::getline(is,line);
			for (size_t j=1;j<=nstates;j++)
			{
				std::getline(is,line);
				std::vector<std::string> tokens = split(line,' ');
				size_t row_idx = std::stoul(tokens[0]);
				for(size_t k=1;k<tokens.size();k++)
				{
					size_t col_idx = k-1+(i-1)*5;
					if(col_idx>=ZERO_ORDER_H.cols() || row_idx-1 >= ZERO_ORDER_H.rows())
						opencap_throw("Error: State index of out bounds. There is a problem with the OpenMolcas output file. Exiting...");
					ZERO_ORDER_H(row_idx-1,col_idx)=std::stod(tokens[k]);
				}
			}
		}
	}
	else
	{
    	opencap_throw("Error: I couldn't read:" + filename);
	}
	return ZERO_ORDER_H;
}


std::vector<Atom> read_geometry_from_rassi(std::string filename)
{
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	auto nsym = file.readAttribute<long>("NSYM", "/");
	std::string atom_nums_tag,coords_tag;
	if (nsym>1)
	{
		atom_nums_tag="DESYM_CENTER_ATNUMS";
		coords_tag = "DESYM_CENTER_COORDINATES";
	}
	else
	{
		atom_nums_tag="CENTER_ATNUMS";
		coords_tag = "CENTER_COORDINATES";
	}
	std::vector<Atom> atoms;
	std::vector<long> atomic_nums;
	Eigen::MatrixXd coords_mat;
	file.readDataset(atomic_nums,atom_nums_tag);
	file.readDataset(coords_mat,coords_tag);
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
		int l = basis_ids(i,2);
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
		int l = prim_ids(i,1);
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
