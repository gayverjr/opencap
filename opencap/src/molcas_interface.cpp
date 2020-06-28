/*
 * readMolcasHDF5.cpp
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */
#include <armadillo>
#include <algorithm>
#include <iostream>
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

std::array<std::vector<std::vector<arma::mat>>,2> read_rassi_tdms(std::string filename)
{
	std::vector<std::vector<arma::mat>> alpha_opdms;
	std::vector<std::vector<arma::mat>> beta_opdms;

	std::vector<std::vector<Eigen::MatrixXd>> alpha_opdms2;
	std::vector<std::vector<Eigen::MatrixXd>> beta_opdms2;

	//first lets load in the hdf5 file
	//arma::cube rassi_data, spin_density;
	Eigen::Tensor<double,3> rass_data, spin_dens;
	h5pp::File file(filename, h5pp::FilePermission::READONLY);
	file.readDataset(rass_data,"SFS_TRANSITION_DENSITIES");
	file.readDataset(spin_dens,"SFS_TRANSITION_SPIN_DENSITIES");
	arma::cube rassi_data, spin_density;
	rassi_data.load(arma::hdf5_name(filename, "SFS_TRANSITION_DENSITIES"));
	spin_density.load(arma::hdf5_name(filename, "SFS_TRANSITION_SPIN_DENSITIES"));
	//now lets loop over the matrices
	for (size_t i=0;i<rassi_data.n_slices;i++)
	{
		std::vector<arma::mat> alpha_state_row;
		std::vector<arma::mat> beta_state_row;
		for(size_t j=0;j<rassi_data.n_slices;j++)
		{
			arma::mat opdm = rassi_data.slice(i).col(j);
			arma::mat spin = spin_density.slice(i).col(i);
			spin.reshape(sqrt(spin.n_rows),sqrt(spin.n_rows));
			opdm.reshape(sqrt(opdm.n_rows),sqrt(opdm.n_rows));
			alpha_state_row.push_back(0.5*(opdm+spin));
			beta_state_row.push_back(0.5*(opdm-spin));
		}
		alpha_opdms.push_back(alpha_state_row);
		beta_opdms.push_back(beta_state_row);
	}
	//symmetric
    for (size_t i=0;i<rassi_data.n_slices;i++)
    {
    	for(size_t j=0;j<i;j++)
    	{
    		alpha_opdms[i][j]= alpha_opdms[j][i];
    		beta_opdms[i][j]= beta_opdms[j][i];
    	}
    }

    const auto& d = rass_data.dimensions();
    cout << "Dim size: " << d.size() << ", dim 0: " << d[0]
         << ", dim 1: " << d[1] << ", dim2:" << d[2] << std::endl;
    for(long i=0;i<d[0];i++)
    {
		std::vector<Eigen::MatrixXd> alpha_state_row;
		std::vector<Eigen::MatrixXd> beta_state_row;
		for (long j=0;j<d[1];j++)
		{
			long nbas = sqrt(d[2]);
			Eigen::array<long,3> offset = {i,j,0};         //Starting point
			Eigen::array<long,3> extent = {i+1,j+1,nbas}; //end point
			Eigen::Tensor<double, 2> dmt_slice = rass_data.slice(offset, extent).reshape(Eigen::array<long,2>{nbas,nbas});
			Eigen::MatrixXd  dmt_mat =  Tensor_to_Matrix(dmt_slice,nbas,nbas);
			Eigen::Tensor<double, 2> spin_slice = spin_dens.slice(offset, extent).reshape(Eigen::array<long,2>{nbas,nbas});
			Eigen::MatrixXd  spin_mat =  Tensor_to_Matrix(dmt_slice,nbas,nbas);
			Eigen::MatrixXd alpha_opdm = 0.5*(dmt_mat+spin_mat);
			Eigen::MatrixXd beta_opdm = 0.5*(dmt_mat-spin_mat);
			alpha_state_row.push_back(alpha_opdm);
			beta_state_row.push_back(beta_opdm);
		}
		alpha_opdms2.push_back(alpha_state_row);
		beta_opdms2.push_back(beta_state_row);
    }
    for (size_t i=0;i<d[0];i++)
    {
    	for(size_t j=0;j<i;j++)
    	{
    		alpha_opdms2[i][j]= alpha_opdms2[j][i];
    		beta_opdms2[i][j]= beta_opdms2[j][i];
    	}
    }

    return {alpha_opdms,beta_opdms};
}

arma::mat read_rassi_overlap(std::string filename)
{
	arma::mat overlap_mat;
	overlap_mat.load(arma::hdf5_name(filename, "AO_OVERLAP_MATRIX"));
	overlap_mat.reshape(sqrt(overlap_mat.n_cols),sqrt(overlap_mat.n_cols));
	return overlap_mat;
}

arma::mat read_mscaspt2_heff(size_t nstates, std::string filename)
{
	arma::mat ZERO_ORDER_H(nstates,nstates);
	ZERO_ORDER_H.zeros();
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
