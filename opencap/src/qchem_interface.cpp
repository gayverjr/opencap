/*Copyright (c) 2020 James Gayvert

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
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "qchem_interface.h"
#include "opencap_exception.h"
#include "gto_ordering.h"
#include "utils.h"
#include "BasisSet.h"
#include <Eigen/Dense>


std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_dms(std::string fchk_filename,BasisSet &bs)
{
	std::ifstream is(fchk_filename);
	size_t nstates = 0;
	bool sep_alpha_beta;
	std::string line, rest;
	while(is.peek() != EOF )
	{
    	std::getline(is, line);
    	if(line.find("Alpha  State Density")!= std::string::npos)
    	{
    		sep_alpha_beta=true;
			int num_elements = stoi(split(line,' ').back());
			if(sqrt(num_elements)!=bs.Nbasis)
				opencap_throw("Error: dimensions of DMs do not match specified basis set.");
    		break;
    	}
    	else if(line.find("State Density")!= std::string::npos)
    	{
    		sep_alpha_beta=false;
			int num_elements = stoi(split(line,' ').back());
			if(sqrt(num_elements)!=bs.Nbasis)
				opencap_throw("Error: dimensions of state DMs do not match specified basis set.");
    		break;
    	}
	}
	is.seekg (0, ios::beg);
	while(is.peek()!=EOF)
	{
    	std::getline(is, line);
    	if(line.find("Alpha  State Density")!= std::string::npos && sep_alpha_beta)
    		nstates++;
    	else if(line.find("State Density")!= std::string::npos && !sep_alpha_beta)
    		nstates++;
	}
	if(nstates==0)
		opencap_throw("Error:Unable to find any densities in:" +fchk_filename);
	if(sep_alpha_beta)
		return qchem_read_in_dms_open_shell(fchk_filename,nstates,bs);
	else
		return qchem_read_in_dms_closed_shell(fchk_filename,nstates,bs);


}

std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_in_dms_open_shell(std::string dmat_filename,
		size_t nstates, BasisSet bs)
{
	std::vector< std::vector<Eigen::MatrixXd>> alpha_opdms(nstates, std::vector<Eigen::MatrixXd> (nstates));
	std::vector< std::vector<Eigen::MatrixXd>> beta_opdms(nstates, std::vector<Eigen::MatrixXd> (nstates));
	//start with state density matrices, alpha and beta densities
	std::ifstream is(dmat_filename);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	for (size_t i=0;i<nstates;i++)
    	{
    			//alpha first, then beta
    			for(size_t spin=1;spin<=2;spin++)
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
					if(sqrt(num_elements)!=bs.Nbasis)
						opencap_throw("Error: dimensions of TDMs do not match specified basis set.");
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
					Eigen::MatrixXd st_opdm(bs.Nbasis,bs.Nbasis);
					st_opdm=Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
					fill_mat(matrix_elements,st_opdm);
					to_opencap_ordering(st_opdm,bs,get_qchem_ids(bs));
					if(spin==1)
						alpha_opdms[i][i]=st_opdm;
					else
						beta_opdms[i][i]=st_opdm;
    			}
		}
    	//now tdms
    	for(size_t i=0;i<nstates;i++)
    	{
    		for (size_t j=i+1;j<nstates;j++)
    		{
    			//alpha first, then beta
    			for (size_t spin=1;spin<=2;spin++)
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
					if(sqrt(num_elements)!=bs.Nbasis)
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
					st_opdm=Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
					fill_mat(matrix_elements,st_opdm);
					to_opencap_ordering(st_opdm,bs,get_qchem_ids(bs));
					if(spin==1)
					{
						alpha_opdms[i][j]=st_opdm;
						alpha_opdms[j][i]=st_opdm;
					}
					else
					{
						beta_opdms[i][j]=st_opdm;
						beta_opdms[j][i]=st_opdm;
					}
    			}
    		}
    	}
    }
    return {alpha_opdms,beta_opdms};
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
		size_t lines_to_read = num_elements%5==0 ? (num_elements/5) : num_elements/5+1;
		std::vector<double> matrix_elements;
		for (size_t k=1;k<=lines_to_read;k++)
		{
			std::getline(is,line);
			std::vector<std::string> tokens = split(line,' ');
			for (auto token:tokens)
				matrix_elements.push_back(std::stod(token));
		}
		fill_LT(matrix_elements,smat);
		to_opencap_ordering(smat,bs,get_qchem_ids(bs));
    }
    return smat;
}

Eigen::MatrixXd read_qchem_energies(size_t nstates,std::string method,std::string output_file)
{
	Eigen::MatrixXd ZERO_ORDER_H(nstates,nstates);
	ZERO_ORDER_H=Eigen::MatrixXd::Zero(nstates,nstates);
	transform(method.begin(),method.end(),method.begin(),::toupper);
	std::ifstream is(output_file);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	size_t state_idx = 1;
    	while (state_idx<=nstates)
    	{
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before "+ std::to_string(nstates) + " energies were found.");
    		std::string line_to_find = method +" transition " + std::to_string(state_idx);
    		if (line.find(line_to_find)!= std::string::npos)
    		{
				std::getline(is,line);
				ZERO_ORDER_H(state_idx-1,state_idx-1) = std::stod(split(line,' ')[3]);
				state_idx++;
    		}
    		else
    			std::getline(is,line);
    	}
    }
    return ZERO_ORDER_H;
}

std::array<std::vector<std::vector<Eigen::MatrixXd>>,2> qchem_read_in_dms_closed_shell(std::string dmat_filename,
		size_t nstates, BasisSet bs)
{
	std::vector< std::vector<Eigen::MatrixXd>> opdms(nstates, std::vector<Eigen::MatrixXd> (nstates));
	//start with state density matrices
	std::ifstream is(dmat_filename);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	for (size_t i=0;i<nstates;i++)
    	{
				while(line.find("State Density")== std::string::npos)
				{
					std::getline(is,line);
    	    		if (is.peek()==EOF)
    	    			opencap_throw("Error: Reached end of file before densities for "+
    	    					std::to_string(nstates) + " states were found.");
				}
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
				Eigen::MatrixXd st_opdm(bs.Nbasis,bs.Nbasis);
				st_opdm=Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
				fill_mat(matrix_elements,st_opdm);
				to_opencap_ordering(st_opdm,bs,get_qchem_ids(bs));
				opdms[i][i]=0.5*st_opdm;
		}
    	//now tdms
    	for(size_t i=0;i<nstates;i++)
    	{
    		for (size_t j=i+1;j<nstates;j++)
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
				st_opdm=Eigen::MatrixXd::Zero(bs.Nbasis,bs.Nbasis);
				fill_mat(matrix_elements,st_opdm);
				to_opencap_ordering(st_opdm,bs,get_qchem_ids(bs));
				opdms[i][j]=0.5*st_opdm;
				opdms[j][i]=0.5*st_opdm;
    		}
    	}
    }
    return {opdms,opdms};

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
    {
    	opencap_throw("Error: I couldn't read:" + fchk_filename);
    }
    return atoms;
}

BasisSet read_basis_from_fchk(std::string fchk_filename, std::vector<Atom> atoms)
{
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
    return bs;
}
