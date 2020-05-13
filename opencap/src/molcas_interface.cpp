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

std::array<std::vector<std::vector<arma::mat>>,2> read_rassi_tdms(std::string filename)
{
	std::vector<std::vector<arma::mat>> alpha_opdms;
	std::vector<std::vector<arma::mat>> beta_opdms;
	//first lets load in the hdf5 file
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
	std::cout << "Successfully read in densities from file:" << filename << std::endl;
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
	std::cout << "Successfully read in zeroth order Hamiltonian from file:" << filename << std::endl;
	return ZERO_ORDER_H;
}
