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
 * qchem_interface.cpp
 */

#include "qchem_interface.h"

#include <algorithm>
#include <chrono>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "gto_ordering.h"
#include "opencap_exception.h"
#include "utils.h"


void qchem_parse_fchk_dms(std::string dmat_filename,std::vector<std::vector<Eigen::MatrixXd>> &alpha_opdms,
		std::vector<std::vector<Eigen::MatrixXd>> &beta_opdms, size_t nstates, size_t ntdm, BasisSet &bs,
		bool symmetric_rdm, bool do_spin)
{
	alpha_opdms = std::vector< std::vector<Eigen::MatrixXd>>(nstates, std::vector<Eigen::MatrixXd> (nstates));
	beta_opdms = std::vector< std::vector<Eigen::MatrixXd>>(nstates, std::vector<Eigen::MatrixXd> (nstates));
	int lt_number = bs.Nbasis*(bs.Nbasis+1)/2;
	int lt_tdm = nstates*(nstates+1)/2 - nstates;
	std::ifstream is(dmat_filename);
	std::vector<Eigen::MatrixXd> alpha_tdms;
	std::vector<Eigen::MatrixXd> beta_tdms;
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	auto start = std::chrono::high_resolution_clock::now();
    	for (size_t i=0;i<nstates;i++)
    	{
    			//alpha first, then beta
    			for(size_t spin=0;spin<=do_spin;spin++)
    			{
					while(line.find("State Density")== std::string::npos)
					{
						std::getline(is,line);
	    	    		if (is.peek()==EOF)
	    	    			opencap_throw("Error: Reached end of file before densities for "+
	    	    					std::to_string(nstates) + " states were found.");
					}
					//last part of line should be number of elements to read
					int num_elements = stoi(split(line,' ').back());
					if(sqrt(num_elements)!=bs.Nbasis && num_elements!=lt_number)
						opencap_throw("Error: dimensions of TDMs do not match specified basis set.");
					size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
					std::vector<double> matrix_elements;
					for (size_t k=1;k<=lines_to_read;k++)
					{
						std::getline(is,line);
						std::vector<std::string> tokens = split(line,' ');
						for (auto token:tokens)
							matrix_elements.push_back(std::stod(token));
					}
					Eigen::MatrixXd st_opdm(bs.Nbasis,bs.Nbasis);
					if(symmetric_rdm)
						fill_LT<double>(matrix_elements,st_opdm);
					else
						fill_mat<double>(matrix_elements,st_opdm);
					to_opencap_ordering(st_opdm,bs,get_qchem_ids(bs));
					if(spin==0)
					{
						alpha_opdms[i][i]=st_opdm;
						if(!do_spin)
							beta_opdms[i][i]=st_opdm;
					}
					else
						beta_opdms[i][i]=st_opdm;
    			}
		}
        is.seekg (0, ios::beg);
    	for(size_t i=0;i<ntdm;i++)
    	{
			//alpha first, then beta
			for (size_t spin=0;spin<=do_spin;spin++)
			{
				while(line.find("Transition DM")== std::string::npos)
				{
					std::getline(is,line);
					if (is.peek()==EOF)
						opencap_throw("Error: Reached end of file before densities for "+
								std::to_string(nstates) + " states were found.");
				}
				//last part of line should be number of elements to read
				size_t num_elements = stoi(split(line,' ').back());
				if(sqrt(num_elements)!=bs.Nbasis && num_elements!=lt_number)
					opencap_throw("Error: dimensions of TDMs do not match specified basis set.");
				size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
				std::vector<double> matrix_elements;
				for (size_t k=1;k<=lines_to_read;k++)
				{
					std::getline(is,line);
					std::vector<std::string> tokens = split(line,' ');
					for (auto token:tokens)
						matrix_elements.push_back(std::stod(token));
				}
				Eigen::MatrixXd st_opdm(bs.Nbasis,bs.Nbasis);
				fill_mat<double>(matrix_elements,st_opdm);
				to_opencap_ordering(st_opdm,bs,get_qchem_ids(bs));
				if(spin==0)
				{
					alpha_tdms.push_back(st_opdm);
					if(!do_spin)
						beta_tdms.push_back(st_opdm);
				}
				else
					beta_tdms.push_back(st_opdm);
			}
    	}
    }
    is.close();
    if(ntdm == lt_tdm)
    {
        std::cout << "Warning: TDM M-->N is assumed to be conjugate transpose of "
        << "TDM N-->M where M>N" << std::endl;
    	size_t dm_idx = 0;
    	for(size_t i=0;i<nstates;i++)
    	{
    		for(size_t j=i+1;j<nstates;j++)
    		{
    			alpha_opdms[i][j] = alpha_tdms[dm_idx];
    			alpha_opdms[j][i] = alpha_tdms[dm_idx].adjoint();
    			beta_opdms[i][j] = beta_tdms[dm_idx].adjoint();
    			beta_opdms[j][i] = beta_tdms[dm_idx];
				dm_idx++;
    		}
    	}
    }
    else
    {
    	size_t dm_idx = 0;
    	for(size_t i=0;i<nstates;i++)
    	{
    		for(size_t j=0;j<nstates;j++)
    		{
    			if(i!=j)
    			{
    				alpha_opdms[i][j] = alpha_tdms[dm_idx];
    				beta_opdms[i][j] = beta_tdms[dm_idx];
    				dm_idx++;
    			}
    		}
    	}
    }
}

void qchem_read_dms(std::vector<std::vector<Eigen::MatrixXd>> &alpha_dms,
		std::vector<std::vector<Eigen::MatrixXd>> &beta_dms,
		std::string fchk_filename, BasisSet &bs,size_t num_states)
{
	std::ifstream is(fchk_filename);
	size_t nstates = 0;
	size_t num_tdm = 0;
	bool sep_alpha_beta;
	bool symmetric_rdm = false;
	bool symmetric_tdm = false;
	std::string line, rest;
	if (is.good())
	{
		while(is.peek() != EOF )
		{
			std::getline(is, line);
			if(line.find("State Density")!= std::string::npos)
			{
				if(line.find("Alpha")!=std::string::npos)
				{
					sep_alpha_beta=true;
					int num_elements = stoi(split(line,' ').back());
					if(sqrt(num_elements)!=bs.Nbasis)
					{
						int lt_number = bs.Nbasis*(bs.Nbasis+1)/2;
						if(lt_number == num_elements)
							symmetric_rdm = true;
						else
							opencap_throw("Error: dimensions of DMs do not match specified basis set.");
					}
					break;
				}
				else
				{
					sep_alpha_beta=false;
					int num_elements = stoi(split(line,' ').back());
					if(sqrt(num_elements)!=bs.Nbasis)
					{
						int lt_number = bs.Nbasis*(bs.Nbasis+1)/2;
						if (lt_number == num_elements)
							symmetric_rdm = true;
						else
							opencap_throw("Error: dimensions of DMs do not match specified basis set.");
					}
					break;
				}
			}
		}
		is.seekg (0, ios::beg);
		while(is.peek()!=EOF)
		{
			std::getline(is, line);
			if((line.find("Alpha")!= std::string::npos && line.find("State Density")!= std::string::npos) && sep_alpha_beta)
				nstates++;
			else if(line.find("State Density")!= std::string::npos && !sep_alpha_beta)
				nstates++;
			if((line.find("Alpha")!= std::string::npos && line.find("Transition DM")!= std::string::npos) && sep_alpha_beta)
				num_tdm++;
			else if(line.find("Transition DM")!= std::string::npos && !sep_alpha_beta)
				num_tdm++;
		}
		is.close();
		if(nstates==0)
		{
			opencap_throw("Error: Unable to find any densities in:" +fchk_filename);
		}
		else if(nstates!=num_states)
		{
			opencap_throw("Error: number of states found: " + std::to_string(nstates) + ", does not match "
					"number of states specified in the input: " + std::to_string(num_states)+ ". Exiting...");
		}
		qchem_parse_fchk_dms(fchk_filename,alpha_dms,beta_dms,nstates,num_tdm,bs,symmetric_rdm,sep_alpha_beta);
	}
	else
	    opencap_throw("Error: I couldn't read:" + fchk_filename);
}

Eigen::MatrixXd qchem_read_overlap(std::string dmat_filename, BasisSet bs)
{
	size_t num_bf = bs.Nbasis;
    std::ifstream is(dmat_filename);
	Eigen::MatrixXd smat(num_bf,num_bf);
	smat=Eigen::MatrixXd::Zero(num_bf,num_bf);
    if (is.good())
    {
    	std::string line, rest;
    	while (line.find("Overlap Matrix")== std::string::npos)
    	{
        	std::getline(is, line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before Overlap Matrix was found.");
    	}
    	size_t num_elements = stoi(split(line,' ').back());
    	int lt_number = bs.Nbasis*(bs.Nbasis+1)/2;
    	if(num_elements!=lt_number)
    		opencap_throw("Error: Dimensions of overlap matrix do not match basis.");
		size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
		std::vector<double> matrix_elements;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				matrix_elements.push_back(std::stod(token));
		}
		fill_LT<double>(matrix_elements,smat);
		to_opencap_ordering(smat,bs,get_qchem_ids(bs));
    }
    else
    	opencap_throw("Error: I couldn't read:" + dmat_filename);
    return smat;
}


Eigen::MatrixXd read_qchem_tddft_energies(size_t nstates,std::string output_file)
{
	Eigen::MatrixXd ZERO_ORDER_H(nstates,nstates);
	ZERO_ORDER_H=Eigen::MatrixXd::Zero(nstates,nstates);
	std::ifstream is(output_file);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	size_t state_idx = 1;
    	while (state_idx<nstates)
    	{
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before "+ std::to_string(nstates) + " energies were found. "
    					"Only " + std::to_string(state_idx) + " states were found. Exiting...");
    		if (line.find("Total energy in the final basis set")!= std::string::npos)
    			ZERO_ORDER_H(0,0) = std::stod(split(line,' ')[8]);
			if (line.find("Total energy for state")!= std::string::npos)
			{
					ZERO_ORDER_H(state_idx,state_idx) = std::stod(split(line,' ')[5]);
					state_idx++;
					std::getline(is,line);
			}
    		else
    			std::getline(is,line);
    	}
    }
    else
    	opencap_throw("Error: I couldn't read:" + output_file);
    return ZERO_ORDER_H;
}

Eigen::MatrixXd read_qchem_eom_energies(size_t nstates,std::string output_file)
{
	Eigen::MatrixXd ZERO_ORDER_H(nstates,nstates);
	ZERO_ORDER_H=Eigen::MatrixXd::Zero(nstates,nstates);
	std::ifstream is(output_file);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	size_t state_idx = 1;
    	while (state_idx<=nstates)
    	{
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before "+ std::to_string(nstates) + " energies were found. "
    					"Only " + std::to_string(state_idx) + " states were found. Exiting...");
			if (line.find("Total energy")!= std::string::npos && line.find("Excitation energy")!= std::string::npos)
			{
					ZERO_ORDER_H(state_idx-1,state_idx-1) = std::stod(split(line,' ')[3]);
					state_idx++;
					std::getline(is,line);
			}
    		else
    			std::getline(is,line);
    	}
    }
    else
    	opencap_throw("Error: I couldn't read:" + output_file);
    return ZERO_ORDER_H;
}

std::vector<Atom> read_geometry_from_fchk(std::string fchk_filename)
{
	std::vector<Atom> atoms;
	std::vector<size_t> atom_nums;
	std::vector<double> coords;
	std::ifstream is(fchk_filename);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
		while(line.find("Atomic numbers")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before atomic numbers were found.");
		}
		size_t num_elements = stoi(split(line,' ').back());
		size_t lines_to_read = num_elements%6==0 ? (num_elements/6) : num_elements/6+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				atom_nums.push_back(std::stoi(token));
		}
		while(line.find("Current cartesian coordinates")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before Current Cartesian Coordinates were found.");
		}
		num_elements = stoi(split(line,' ').back());
		lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				coords.push_back(std::stod(token));
		}
		//ok now lets populate our atoms
		for(size_t i=0;i<atom_nums.size();i++)
			atoms.push_back(Atom(atom_nums[i],coords[i*3],coords[i*3+1],coords[i*3+2]));
    }
    else
    	opencap_throw("Error: I couldn't read:" + fchk_filename);
    return atoms;
}

BasisSet read_basis_from_fchk(std::string fchk_filename, std::vector<Atom> atoms)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<int> shell_types;
	std::vector<int> prims_per_shell;
	std::vector<int> atom_ids;
	std::vector<double> exps;
	std::vector<double> coeffs;
	// needed for SP specification
	std::vector<double> p_coeffs;
	bool SP_basis_function = true;
	BasisSet bs;
	for(auto atm:atoms)
		bs.centers.push_back(atm.coords);
	std::vector<shell_id> ids;
	std::ifstream is(fchk_filename);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	//First lets figure out if there are SP functions
    	//shell types
		while(line.find("Shell types")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before shell types were found.");
		}
		size_t num_elements = stoi(split(line,' ').back());
		size_t lines_to_read = num_elements%6==0 ? (num_elements/6) : num_elements/6+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
			{
				shell_types.push_back(std::stoi(token));
			}
		}
		//prims per shell
		while(line.find("Number of primitives per shell")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before Number of primitives per shell were found.");
		}
		num_elements = stoi(split(line,' ').back());
		lines_to_read = num_elements%6==0 ? (num_elements/6) : num_elements/6+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				prims_per_shell.push_back(std::stoi(token));
		}
		//atom ids
		while(line.find("Shell to atom map")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before Shell to atom map was found.");
		}
		num_elements = stoi(split(line,' ').back());
		lines_to_read = num_elements%6==0 ? (num_elements/6) : num_elements/6+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				atom_ids.push_back(std::stoi(token));
		}
		//prims
		while(line.find("Primitive exponents")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before Primitive exponents were found.");
		}
		num_elements = stoi(split(line,' ').back());
		lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				exps.push_back(std::stod(token));
		}
		//coeffs
		while(line.find("Contraction coefficients")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before Contraction coefficients were found.");
		}
		num_elements = stoi(split(line,' ').back());
		lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				coeffs.push_back(std::stod(token));
		}
		while(line.find("P(S=P) Contraction coefficients")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    		{
    			SP_basis_function = false;
    			break;
    		}
		}
		if (SP_basis_function)
		{
			num_elements = stoi(split(line,' ').back());
			lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
			for (size_t k=1;k<=lines_to_read;k++)
			{
				std::getline(is,line);
				std::vector<std::string> tokens = split(line,' ');
				for (auto token:tokens)
					p_coeffs.push_back(std::stod(token));
			}
		}
    }
    else
    	opencap_throw("Error: I couldn't read:" + fchk_filename);
    size_t prim_idx=0;
    for(size_t i=0;i<shell_types.size();i++)
    {
    	//SP
    	if(shell_types[i]==-1)
    	{
    		if(p_coeffs.size()==0)
    			opencap_throw("Error: missing section P(S=P) Contraction coefficients.");
    		Shell s_shell(0,atoms[atom_ids[i]-1].coords);
    		Shell p_shell(1,atoms[atom_ids[i]-1].coords);
			int num_prims = prims_per_shell[i];
			for(int j=1;j<=num_prims;j++)
			{
				s_shell.add_primitive(exps[prim_idx],coeffs[prim_idx]);
				p_shell.add_primitive(exps[prim_idx],p_coeffs[prim_idx]);
				prim_idx++;
			}
			bs.add_shell(s_shell);
			bs.add_shell(p_shell);
    	}
    	else
    	{
			Shell new_shell(abs(shell_types[i]),atoms[atom_ids[i]-1].coords);
			if(shell_types[i]>0 && new_shell.l>1)
				new_shell.pure=false;
			int num_prims = prims_per_shell[i];
			for(int j=1;j<=num_prims;j++)
			{
				new_shell.add_primitive(exps[prim_idx],coeffs[prim_idx]);
				prim_idx++;
			}
			bs.add_shell(new_shell);
    	}
    }
    bs.normalize();
	auto stop = std::chrono::high_resolution_clock::now();
	auto total_time = std::chrono::duration<double>(stop-start).count();
    return bs;
}
