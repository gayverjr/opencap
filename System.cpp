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
#include "molcas_transform.h"
#include "readMolcasHDF5.h"
#include "CAP.h"
#include "overlap.h"
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
			EOMCAP(row_idx,col_idx) =  arma::trace(alpha_dms[row_idx][col_idx]*AO_CAP_MAT)+
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
	else if(pkg=="molcas")
		to_molcas_ordering(AO_CAP_MAT,bs,atoms);
	else
		std::cout << "Other things coming soon!" << std::endl;
}

void System::compute_cap_matrix()
{
	std::cout << "Calculating CAP matrix in AO basis..." << std::endl;
	CAP cap_integrator(atoms,parameters);
	arma::mat cap_mat(bs.num_carts(),bs.num_carts());
	cap_mat.zeros();
	auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_cap_mat(cap_mat,bs);
	auto stop = std::chrono::high_resolution_clock::now();
	std::cout << "Integration time:" << std::chrono::duration<double, std::milli>(stop-start).count() << std::endl;
	uniform_cart_norm(cap_mat,bs);
	arma::mat cap_spherical(bs.Nbasis,bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,bs);
	AO_CAP_MAT = cap_spherical;
	//re-order cap matrix based on electronic structure package dms came from
	reorder_cap();
	compute_cap_correlated_basis();

}

void System::ang_to_bohr()
{
	for (auto atm:atoms)
		atm.ang_to_bohr();
}

bool System::check_overlap_matrix()
{
	//get overlap matrix
	arma::mat Smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,Smat);
	uniform_cart_norm(Smat,bs);
	arma::mat spherical_ints(bs.Nbasis,bs.Nbasis);
	cart2spherical(Smat,spherical_ints,bs);
	if (parameters["package"]=="qchem")
	{
		to_molden_ordering(spherical_ints, bs);
		auto qchem_smat = qchem_read_overlap(parameters["fchk_file"],bs.Nbasis);
		bool conflicts = false;
		for (size_t i=0;i<qchem_smat.n_rows;i++)
		{
			for(size_t j=0;j<qchem_smat.n_cols;j++)
			{
				if (abs(qchem_smat(i,j)-spherical_ints(i,j))>1E-5)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
					std::cout << "qchem says:" << qchem_smat(i,j) << std::endl;
					std::cout << "OpenCAP says:" << spherical_ints(i,j) << std::endl;
					std::cout << abs(qchem_smat(i,j)-spherical_ints(i,j)) << std::endl;
					conflicts = true;
				}
			}
		}
		if (conflicts)
		{
			std::cout << "Should probably throw an error here." << std::endl;
			return false;
		}
		else
			return true;
	}

	else if (parameters["package"]=="molcas")
	{
		to_molcas_ordering(spherical_ints,bs,atoms);
		arma::mat overlap_mat;
		overlap_mat.load(arma::hdf5_name(parameters["rassi_h5"], "AO_OVERLAP_MATRIX"));
		overlap_mat.reshape(sqrt(overlap_mat.n_cols),sqrt(overlap_mat.n_cols));
	    std::cout << std::fixed << std::setprecision(10);
		bool conflicts = false;
		for (size_t i=0;i<overlap_mat.n_rows;i++)
		{
			for(size_t j=0;j<overlap_mat.n_cols;j++)
			{
				if (abs(overlap_mat(i,j)-spherical_ints(i,j))>1E-12)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
					std::cout << "Molcas says:" << overlap_mat(i,j) << std::endl;
					std::cout << "OpenCAP says:" << spherical_ints(i,j) << std::endl;
					std::cout << abs(overlap_mat(i,j)-spherical_ints(i,j)) << std::endl;
					conflicts = true;
				}
			}
		}
		if (!conflicts)
			return true;
		else
			return false;
	}
	else
	{
		std::cout << "WARNING: Unable to check overlap matrix, but continuing anyways..." << std::endl;
		return true;
	}
}

bool System::read_in_dms()
{

	stringstream ss(parameters["nstates"]);
	size_t nstates;
	ss >> nstates;
	if (parameters["package"]=="qchem")
	{
		std::string dmat_filename = parameters["fchk_file"];
		auto parsed_dms = qchem_read_in_dms(dmat_filename,nstates,bs.Nbasis);
		alpha_dms = parsed_dms[0];
		beta_dms = parsed_dms[1];
		return true;
	}
	else if (parameters["package"]=="molcas")
	{
		auto parsed_dms = read_rassi_HDF5(parameters["rassi_h5"]);
		alpha_dms = parsed_dms[0];
		beta_dms = parsed_dms[1];
		return true;
	}
	else
		std::cout << "Only q-chem and molcas are supported." << std::endl;
	return false;

}

bool System::read_molcas_Heff()
{
	return true;
}

bool System::read_qchem_energies()
{
	ZERO_ORDER_H = arma::mat(nstates,nstates);
	ZERO_ORDER_H.zeros();
	std::string method = parameters["method"];
	transform(method.begin(),method.end(),method.begin(),::toupper);
	std::cout << "Reading file:" << parameters["qc_output"] << std::endl;
	std::ifstream is(parameters["qc_output"]);
    if (is.good())
    {
    	std::string line, rest;
    	std::getline(is, line);
    	size_t state_idx = 1;
    	while (state_idx<=nstates)
    	{
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
    std::cout << "Printing Zeroth Order Hamiltonian" << std::endl;
    ZERO_ORDER_H.raw_print();
    return false;
}

bool System::read_in_zero_order_H()
{
	if (parameters["package"]=="qchem")
		return read_qchem_energies();
	else if (parameters["package"]=="molcas")
		return read_molcas_Heff();
	else
		std::cout << "Only q-chem and molcas formats are supported." << std::endl;
	return false;

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
	if (read_in_zero_order_H())
		std::cout << "Successfully read in zeroth order Hamiltonian." << std::endl;
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

bool System::verify_method(std::string key)
{
	std::string package_name = parameters[key];
	if(parameters.find("method")==parameters.end())
	{
		std::cout << "Error: missing the 'method' keyword. Please choose a supported package/method." << std::endl;
		return false;
	}
	std::string method = parameters["method"];
	if (package_name=="qchem"||package_name=="q-chem")
	{
		std::vector<std::string> supported = {"eomea","eomee","eomip"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
		{
			std::cout << "Error: unsupported Q-Chem method. OpenCAP currently supports: 'eom-ea-ccsd','eom-ee-ccsd','eom-ip-ccsd'."<< std::endl;
			return false;
		}
		if (parameters.find("fchk_file")==parameters.end())
		{
			std::cout << "Error: missing keyword: fchk_file." << std::endl;
			return false;
		}
		if (parameters.find("qc_output")==parameters.end())
		{
			std::cout << "Error: missing keyword: qc_output" << std::endl;
			return false;
		}
		return true;
	}
	else if(package_name=="molcas")
	{
		std::vector<std::string> supported = {"ms-caspt2","xms-caspt2","qd-nevpt2"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
		{
			std::cout << "Error: unsupported Molcas method. OpenCAP currently supports: 'dmrgscf','nevpt2','caspt2','rasscf'."<< std::endl;
			return false;
		}
		if (parameters.find("rassi_h5")==parameters.end())
		{
			std::cout << "Error: missing keyword: rassi_h5." << std::endl;
			return false;
		}
		if (parameters.find("molcas_output")==parameters.end())
		{
			std::cout << "Error: missing keyword: molcas_output" << std::endl;
			return false;
		}
		return true;
	}
	else
	{
		std::cout << "Error: unsupported package. Only Q-Chem and Molcas are currently supported." << std::endl;
		return false;
	}
}

//check that requirements have been specified, and set optional params to their defaults
bool System::check_parameters()
{
	std::map<std::string, std::string> defaults =
	{{"bohr_coordinates", "false"}, {"radial_precision", "14"}, {"angular_points", "590"},{"cart_bf",""},{"basis_file",""}};
	std::vector<std::string> requirements = {"package","nstates","cap_type","basis_set"};
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
			else if(key=="package")
			{
				if(!verify_method(key))
					return false;
			}
		}
		else
		{
			std::cout << "Error: missing required key \"" << key << "\""<< std::endl;
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
