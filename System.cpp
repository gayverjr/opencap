/*
 * System.cpp
 *
 *  Created on: Feb 18, 2020
 *      Author: JG
 */
#include "System.h"
#include "BasisSet.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "CAP.h"
#include <string>
#include <map>
#include <list>
#include <armadillo>
#include "Atom.h"
#include "read_qchem_fchk.h"
#include "numerical.h"
#include "utils.h"
#include "transforms.h"
#include "molden_transform.h"
#include "CAP.h"
#include <cmath>
#include <limits>


void System::compute_cap_correlated_basis()
{
	arma::mat EOMCAP(nstates,nstates);
	EOMCAP.zeros();
    std::cout << std::fixed << std::setprecision(10);
	for (size_t row_idx=0;row_idx<EOMCAP.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<EOMCAP.n_cols;col_idx++)
		{
			EOMCAP(row_idx,col_idx) = arma::trace(alpha_dms[row_idx][col_idx]*AO_CAP_MAT)+
									  arma::trace(beta_dms[row_idx][col_idx]*AO_CAP_MAT);
			EOMCAP(row_idx,col_idx) = -1.0* EOMCAP(row_idx,col_idx);
		}
	}
	CORRELATED_CAP_MAT = EOMCAP;
}

void System::reorder_cap()
{
	std::string pkg = parameters["package"];
	if (pkg=="q-chem"||pkg=="qchem")
		to_molden_ordering(AO_CAP_MAT, bs);
	else
	{
		std::cout << "NYI" << std::endl;
	}
}

void System::compute_cap_matrix()
{
	std::cout << "Calculating CAP matrix in AO basis...";
	CAP cap_integrator(atoms,parameters);
	arma::mat cap_mat(bs.num_carts(),bs.num_carts());
	cap_mat.zeros();
	auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_cap_mat(cap_mat,bs);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	uniform_cart_norm(cap_mat,bs);
	arma::mat cap_spherical(bs.Nbasis,bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,bs);
	AO_CAP_MAT = cap_spherical;
	//re-order cap matrix based on electronic structure package dms came from
	reorder_cap();
	std::cout << "done." << std::endl;
	compute_cap_correlated_basis();

}

void System::ang_to_bohr()
{
	for (auto atm:atoms)
		atm.ang_to_bohr();
}

bool System::check_overlap_matrix()
{return true;}

bool System::read_in_dms()
{
	stringstream ss(parameters["nstates"]);
	size_t nstates;
	ss >> nstates;
	std::string dmat_filename = parameters["data_file"];
	auto parsed_dms = qchem_read_in_dms(dmat_filename,nstates,bs.Nbasis);
	alpha_dms = parsed_dms[0];
	beta_dms = parsed_dms[1];
	return true;
}

System::System(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	parameters=params;
	std::cout << "Checking input file for required keywords...";
	if(check_parameters())
		std::cout << "ok." << std::endl;
	else
		return;
	atoms = geometry;
	if(parameters["bohr_coordinates"]=="false")
		ang_to_bohr();
	bs = BasisSet(atoms,parameters);
	std::cout << "Checking overlap matrix...";
	if (check_overlap_matrix())
		std::cout << "ok." << std::endl;
	stringstream ss(parameters["nstates"]);
	ss >> nstates;
	std::cout << "Reading in the density matrices... for " <<nstates << " states." << std::endl;
	if (read_in_dms())
		std::cout << "all set!" << std::endl;
}

bool System::verify_cap_parameters(std::string key)
{
	if(parameters[key]=="box" || parameters[key]=="1")
	{
		if(parameters.find("cap_x")==parameters.end())
			return false;
		if(parameters.find("cap_y")==parameters.end())
			return false;
		if (parameters.find("cap_z")==parameters.end())
			return false;
	}
	else if (key=="voronoi" || key=="2")
	{
		if(parameters.find("r_cut")==parameters.end())
			return false;
	}
	else
		return false;
	return true;
}

//check that requirements have been specified, and set optional params to their defaults
bool System::check_parameters()
{
	std::map<std::string, std::string> defaults =
	{{"bohr_coordinates", "false"}, {"radial_precision", "14"}, {"angular_points", "590"},{"cart_bf",""}};
	std::vector<std::string> requirements = {"method","package","data_file","nstates","cap_type","basis_set"};
	//first lets verify that all of our requirements are there
	for (std::string key:requirements)
	{
		if (parameters.find(key)!=parameters.end())
		{
			if (key=="cap_type")
			{
				if(!verify_cap_parameters(key))
					return false;
			}
			if (key=="basis_set" && parameters[key]=="gen")
			{
				if (parameters.find("basis_file")==parameters.end())
				{
					std::cout << "Need to specify a basis file when using gen basis." << std::endl;
					return false;
				}
			}
		}
		else
		{
			std::cout << "Error: missing key \"" << key << "\""<< std::endl;
			return false;
		}
	}
	//now let's set optional parameters to their defaults if not specified
	for (const auto &pair:defaults)
	{
		if(parameters.find(pair.first)==parameters.end())
			parameters[pair.first]=pair.second;
	}
	return true;
}

/*
std::vector<Atom> System::read_xyz(std::string xyz_name)
{
	ifstream xyzfile;
	xyzfile.open(xyz_name);
	//read in number of atoms
	size_t natom;
	xyzfile >> natom;
	std::string rest_of_line;
	std::getline(xyzfile, rest_of_line);
	//read in comment
	std::string comment;
	std::getline(xyzfile, comment);
	//now the coordinates
	std::vector<Atom> atoms;
	for (size_t i=0;i<natom;i++)
	{
	    // read line
	    std::string linestr;
	    std::getline(xyzfile, linestr);
	    std::istringstream iss(linestr);
	    std::string element_symbol;
	    double x, y, z;
	    iss >> element_symbol >> x >> y >> z;
	    transform(element_symbol.begin(),element_symbol.end(),element_symbol.begin(),::toupper);
	    atoms.push_back(Atom(element_symbol,x,y,z));
	}
	return atoms;
}
*/
