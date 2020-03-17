/*
 * read_qchem_fchk.cpp
 *
 *  Created on: Mar 12, 2020
 *      Author: JG
 */
#include <armadillo>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "read_qchem_fchk.h"

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
	  if (!token.empty())
		  tokens.push_back(token);
   }
   return tokens;
}

void fill_mat(std::vector<double> &matrix_elements, arma::mat &opdm)
{
	size_t vec_idx = 0;
	for (size_t row_idx=0;row_idx<opdm.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<opdm.n_cols;col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
	}
}

void fill_LT(std::vector<double> matrix_elements, arma::mat &opdm)
{
	size_t vec_idx = 0;
	size_t row_idx = 0;
	while(row_idx<opdm.n_rows && vec_idx<matrix_elements.size())
	{
		//elements are added to each column <= row index
		for (size_t col_idx=0;col_idx<=row_idx;col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			opdm(col_idx,row_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
		row_idx++;
	}
}

size_t total_TDMs_to_read(size_t nstates)
{
	size_t numTDMs = 0;
	for (size_t i=1;i<nstates;i++)
	{
		for (size_t j=i;j<nstates;j++)
			numTDMs++;
	}
	return numTDMs;
}

size_t get_TDM_start(size_t nstates, size_t state_idx)
{
	size_t tdm_idx = 2*nstates;
	for (size_t i=1;i<state_idx;i++)
	{
		for (size_t j=i;j<nstates;j++)
			tdm_idx+=2;
	}
	return tdm_idx;
}

//currently this is written for open shell systems which have alpha and beta densities
std::array<std::vector<std::vector<arma::mat>>,2> qchem_read_in_dms(std::string dmat_filename,size_t nstates, size_t num_bf)
{
	std::vector<arma::mat> opdms;
	//start with state density matrices, alpha and beta densities
    std::ifstream is(dmat_filename);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	while (opdms.size() < nstates*2)
    	{
				if (line.find("State Density")!= std::string::npos)
				{
					//last part of line should be number of elements to read
					size_t num_elements = stoi(split(line,' ').back());
					size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
					std::vector<double> matrix_elements;
					for (size_t k=1;k<=lines_to_read;k++)
					{
						std::getline(is,line);
						std::vector<std::string> tokens = split(line,' ');
						for (auto token:tokens)
						{
							matrix_elements.push_back(std::stod(token));
						}
					}
					arma::mat st_opdm(num_bf,num_bf);
					st_opdm.zeros();
					fill_LT(matrix_elements,st_opdm);
					opdms.push_back(st_opdm);
				}
				else
					std::getline(is,line);
    	}
    	//now tdms for alpha and beta densities
    	while (opdms.size() < nstates*2 +total_TDMs_to_read(nstates)*2)
    	{
			if (line.find("Transition DM")!= std::string::npos)
			{
				//last part of line should be number of elements to read
				size_t num_elements = stoi(split(line,' ').back());
				size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
				std::vector<double> matrix_elements;
				for (size_t k=1;k<=lines_to_read;k++)
				{
					std::getline(is,line);
					std::vector<std::string> tokens = split(line,' ');
					for (auto token:tokens)
					{
						matrix_elements.push_back(std::stod(token));
					}
				}
				arma::mat st_opdm(num_bf,num_bf);
				st_opdm.zeros();
				fill_mat(matrix_elements,st_opdm);
				opdms.push_back(st_opdm);
			}
			else
				std::getline(is,line);
    	}
    }
    //now that we have our density matrices, lets organize them into a
    //handy form which corresponds to how they'll actually be used
    std::vector<std::vector<arma::mat>> alpha_densities;
    std::vector<std::vector<arma::mat>> beta_densities;
    for (size_t state_idx=1; state_idx<=nstates;state_idx++)
    {
    	std::vector<arma::mat> alpha_row;
    	std::vector<arma::mat> beta_row;
    	//fill with placeholders
    	for (size_t i=1;i<state_idx;i++)
    	{
    		alpha_row.push_back(arma::mat(1,1));
    		beta_row.push_back(arma::mat(1,1));
    	}
    	//state densities
    	alpha_row.push_back(opdms[2*state_idx-2]);
    	beta_row.push_back(opdms[2*state_idx-1]);
    	//transition densities
    	size_t tdm_start = get_TDM_start(nstates,state_idx);
    	size_t tdm_end = tdm_start + 2*(nstates-state_idx);
    	for (size_t j=tdm_start;j<tdm_end;j+=2)
    	{
        	alpha_row.push_back(opdms[j]);
        	beta_row.push_back(opdms[j+1]);
    	}
    	alpha_densities.push_back(alpha_row);
    	beta_densities.push_back(beta_row);
    }
    //fill in the placeholders
    for (size_t i=0;i<nstates;i++)
    {
    	for(size_t j=0;j<i;j++)
    	{
    		alpha_densities[i][j]= alpha_densities[j][i];
    		beta_densities[i][j]= beta_densities[j][i];
    	}
    }
    return {alpha_densities,beta_densities};

}

arma::mat qchem_read_overlap(std::string dmat_filename, size_t num_bf)
{
    std::ifstream is(dmat_filename);
	arma::mat smat;
	smat.zeros(num_bf,num_bf);
    if (is.good())
    {
    	std::string line, rest;
    	while (line.find("Overlap Matrix")== std::string::npos)
        	std::getline(is, line);
    	size_t num_elements = stoi(split(line,' ').back());
		size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
		std::vector<double> matrix_elements;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
			{
				matrix_elements.push_back(std::stod(token));
			}
		}
		fill_LT(matrix_elements,smat);
    }
    return smat;
}



