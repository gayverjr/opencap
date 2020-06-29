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
#include <stdio.h>
#include <h5pp/h5pp.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar,int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}

std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> read_rassi_tdms(std::string filename)
{
	std::vector<std::vector<Eigen::MatrixXd>> alpha_opdms;
	std::vector<std::vector<Eigen::MatrixXd>> beta_opdms;

	//first lets load in the hdf5 file
	Eigen::Tensor<double,3> rass_data, spin_dens;
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
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

Eigen::MatrixXd read_rassi_overlap(std::string filename)
{
	Eigen::VectorXd overlap_mat;
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	file.readDataset(overlap_mat,"AO_OVERLAP_MATRIX");
	int nbas = sqrt(overlap_mat.size());
	Eigen::MatrixXd smat;
	smat = reshape(overlap_mat,nbas,nbas);
	return smat;
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
