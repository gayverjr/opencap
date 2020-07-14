/*
 * System.cpp
 *
 *  Created on: Feb 18, 2020
 *      Author: JG
 */
#include <pybind11/pybind11.h>
#include "System.h"
#include "BasisSet.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <map>
#include <list>
#include "Atom.h"
#include "utils.h"
#include "transforms.h"
#include "gto_ordering.h"
#include "keywords.h"
#include "CAP.h"
#include "overlap.h"
#include <cmath>
#include <limits>
#include "opencap_exception.h"
#include "molcas_interface.h"
#include "qchem_interface.h"
#include "molden_parser.h"
#include <Eigen/Dense>

System::System(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	try
	{
		python = false;
		parameters=params;
		atoms = geometry;
		//set defaults for cart_bf and bohr_coordinates
		std::map<std::string, std::string> defaults = {{"cart_bf", ""}, {"bohr_coordinates", "true"}};
		for (const auto &pair:defaults)
		{
			if(parameters.find(pair.first)==parameters.end())
				parameters[pair.first]=pair.second;
		}
		if(parameters.find("bohr_coordinates")!=parameters.end() && parameters["bohr_coordinates"]=="false")
		{
			for (size_t i=0;i<atoms.size();i++)
				atoms[i].ang_to_bohr();
		}
		bs = BasisSet(atoms,parameters);
		verify_system();
		std::cout << "Number of basis functions:" << bs.Nbasis << std::endl;
		//now construct overlap matrix
		Eigen::MatrixXd Smat(bs.num_carts(),bs.num_carts());
		compute_analytical_overlap(bs,Smat);
		uniform_cart_norm(Smat,bs);
		Eigen::MatrixXd spherical_ints(bs.Nbasis,bs.Nbasis);
		cart2spherical(Smat,spherical_ints,bs);
		OVERLAP_MAT = spherical_ints;
	}
	catch (exception &e)
	{
		opencap_rethrow("Failed to construct System.");
	}

}

System::System(py::dict dict)
{
	python = true;
    for (auto item : dict)
    {
    	std::string key = py::str(item.first).cast<std::string>();
    	std::string value = py::str(item.second).cast<std::string>();
		transform(key.begin(),key.end(),key.begin(),::tolower);
		transform(value.begin(),value.end(),value.begin(),::tolower);
		if(check_keyword(key,"system")||key=="geometry")
			parameters[key]=value;
		else
			opencap_throw("Invalid key:" +key);
    }
    if (parameters.find("molecule")==parameters.end()||parameters.find("basis_file")==parameters.end())
    	opencap_throw("Error: molecule and basis_file keywords are required.");
	if(parameters["molecule"]=="qchem_fchk")
	{
		atoms = read_geometry_from_fchk(parameters["basis_file"]);
		bs = read_basis_from_fchk(parameters["basis_file"],atoms);
	}
	else if (parameters["molecule"]=="molcas_rassi")
	{
		atoms = read_geometry_from_rassi(parameters["basis_file"]);
		bs = read_basis_from_rassi(parameters["basis_file"],atoms);
	}
	else if(parameters["molecule"]=="molden")
	{
		atoms = read_geometry_from_molden(parameters["basis_file"]);
		bs = read_basis_from_molden(parameters["basis_file"],atoms);
	}
	else if(parameters["molecule"]=="read")
	{
		if(parameters.find("geometry")==parameters.end())
			opencap_throw("Error: Need to specify geometry string when molecule is set to \"read.\"");
		atoms = parse_geometry_string(parameters["geometry"]);
		if(parameters.find("bohr_coordinates")!=parameters.end() && parameters["bohr_coordinates"]=="false")
		{
			for (size_t i=0;i<atoms.size();i++)
				atoms[i].ang_to_bohr();
		}
		bs = BasisSet(atoms,parameters);
	}
    verify_system();
	py::print("Number of basis functions:"+std::to_string(bs.Nbasis));
	//now construct overlap matrix
	Eigen::MatrixXd Smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,Smat);
	uniform_cart_norm(Smat,bs);
	Eigen::MatrixXd spherical_ints(bs.Nbasis,bs.Nbasis);
	cart2spherical(Smat,spherical_ints,bs);
	OVERLAP_MAT = spherical_ints;
}

System::System(std::string filename,std::string file_type)
{
	parameters["molecule"]   = filename;
	parameters["basis_file"] = file_type;
	python = false;
	if(file_type=="qchem_fchk")
	{
		atoms = read_geometry_from_fchk(filename);
		bs = read_basis_from_fchk(filename,atoms);
	}
	else if (file_type=="molcas_rassi")
	{
		atoms = read_geometry_from_rassi(filename);
		bs = read_basis_from_rassi(filename,atoms);
	}
	else if(file_type=="molden")
	{
		atoms = read_geometry_from_molden(filename);
		bs = read_basis_from_molden(filename,atoms);
	}
	else
		opencap_throw("That file type isn't supported, sorry.");
	verify_system();
	std::cout << "Number of basis functions:" << bs.Nbasis << std::endl;
	Eigen::MatrixXd Smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,Smat);
	uniform_cart_norm(Smat,bs);
	Eigen::MatrixXd spherical_ints(bs.Nbasis,bs.Nbasis);
	cart2spherical(Smat,spherical_ints,bs);
	OVERLAP_MAT = spherical_ints;
}

std::vector<Atom> System::parse_geometry_string(std::string geometry_string)
{
	std::vector<Atom> atom_list;
	try
	{
		stringstream ss(geometry_string);
		string line;
		while(std::getline(ss,line))
		{
			std::istringstream iss(line);
			std::string element_symbol;
			double x, y, z;
			iss >> element_symbol >> x >> y >> z;
			transform(element_symbol.begin(),element_symbol.end(),element_symbol.begin(),::tolower);
			atom_list.push_back(Atom(element_symbol,x,y,z));
		}
	}
	catch (exception &e)
	{
		opencap_rethrow("Failed to parse geometry string.");
	}
	return atom_list;
}

void System::verify_system()
{
	if(atoms.size()==0)
		opencap_throw("Error: No atoms specified.");
	if(bs.Nbasis==0)
		opencap_throw("Error: Basis set has 0 basis functions.");
}

Eigen::MatrixXd System::get_overlap_mat()
{
	return OVERLAP_MAT;
}

void System::check_overlap_mat(Eigen::MatrixXd smat, std::string ordering, std::string basis_file)
{
	std::vector<bf_id> ids;
	if(ordering=="pyscf")
		ids = get_pyscf_ids(bs);
	else if(ordering=="openmolcas")
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(bs,basis_file);
	}
	else if(ordering=="qchem"||ordering=="molden")
		ids = get_molden_ids(bs);
	else if(ordering=="opencap")
		ids = bs.bf_ids;
	else
		opencap_throw(ordering +" ordering is not supported.");
	to_opencap_ordering(smat,bs,ids);
	if(OVERLAP_MAT.rows() != smat.rows() || OVERLAP_MAT.cols() != smat.cols())
		opencap_throw("Error: Dimension of overlap matrix is incorrect.");
	bool conflicts = false;
	for (size_t i=0;i<smat.rows();i++)
	{
		for(size_t j=0;j<smat.cols();j++)
		{
			if (abs(smat(i,j)-OVERLAP_MAT(i,j))>1E-5)
			{
				std::cout << "Conflict at:" << i << "," << j << std::endl;
				std::cout << ordering + " says:" << smat(i,j) << std::endl;
				std::cout << "OpenCAP says:" << OVERLAP_MAT(i,j) << std::endl;
				conflicts = true;
			}
		}
	}
	if (conflicts)
		opencap_throw("Error: The dimensions of the overlap matrices match, but the elements do not.");
	if(python)
		py::print("Verified overlap matrix.");
	else
		std::cout << "Verified overlap matrix.";
}


