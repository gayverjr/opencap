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
#include <Eigen/Dense>

System::System(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	try
	{
		verify_system_parameters(params);
		parameters=params;
		atoms = geometry;
		if(!bohr_coords())
		{
			for (size_t i=0;i<atoms.size();i++)
				atoms[i].ang_to_bohr();
		}
		bs = BasisSet(atoms,parameters);
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
	std::map<std::string,std::string> params;
    for (auto item : dict)
    {
    	std::string key = py::str(item.first).cast<std::string>();
    	std::string value = py::str(item.second).cast<std::string>();
		transform(key.begin(),key.end(),key.begin(),::tolower);
		transform(value.begin(),value.end(),value.begin(),::tolower);
    	if (key!="geometry")
    		params[key]=value;
    	else
    		set_geometry(value);
    }
    verify_system_parameters(params);
    parameters = params;
	if(!bohr_coords())
	{
		for (size_t i=0;i<atoms.size();i++)
			atoms[i].ang_to_bohr();
	}
	bs = BasisSet(atoms,parameters);
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
	if(file_type=="qchem_fchk")
	{
		atoms = read_atoms_from_fchk(filename);
		bs = read_basis_from_fchk(filename,atoms);
	}
	else if (file_type=="molcas_rassi")
	{
		atoms = read_geometry_from_rassi(filename);
		bs = read_basis_from_rassi(filename,atoms);
	}
	else
		opencap_throw("That file type isn't supported, sorry.");
	Eigen::MatrixXd Smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,Smat);
	uniform_cart_norm(Smat,bs);
	Eigen::MatrixXd spherical_ints(bs.Nbasis,bs.Nbasis);
	cart2spherical(Smat,spherical_ints,bs);
	OVERLAP_MAT = spherical_ints;
	std::cout << "Number of basis functions:" << bs.Nbasis << std::endl;
}

void System::set_geometry(std::string geometry_string)
{
	std::vector<Atom> atom_list;
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
	atoms = atom_list;
}

void System::verify_system_parameters(std::map<std::string, std::string> &params)
{
	if(params.find("basis_file")==params.end())
		opencap_throw("Error: Need to specify a basis set file using the basis_file keyword.");
	std::map<std::string, std::string> defaults = {{"cart_bf", ""}, {"bohr_coordinates", "true"}};
	for (const auto &pair:defaults)
	{
		if(parameters.find(pair.first)==parameters.end())
			parameters[pair.first]=pair.second;
	}
}

bool System::bohr_coords()
{
	std::string my_str = parameters["bohr_coordinates"];
	if (my_str == "false")
		return false;
	else if (my_str == "true")
		return true;
	else if (my_str == "0")
		return false;
	else if (my_str == "1")
		return true;
	else
		opencap_throw("Invalid value for keyword 'bohr_coordinates'");
	return false;
}

Eigen::MatrixXd System::get_overlap_mat(std::string gto_ordering)
{
	Eigen::MatrixXd overlap_copy;
	overlap_copy = OVERLAP_MAT;
	transform(gto_ordering.begin(),gto_ordering.end(),gto_ordering.begin(),::tolower);
	if(gto_ordering=="default")
		return overlap_copy;
	else if (gto_ordering=="pyscf")
	{
		return overlap_copy;
		//return cEigen::MatrixXd_to_arr(overlap_copy);
	}
	else if(gto_ordering=="molden"||gto_ordering=="qchem"||gto_ordering=="q-chem")
	{
		return overlap_copy;
		//return cEigen::MatrixXd_to_arr(overlap_copy);
	}
	else
		opencap_throw("Error: Re-ordering to "+ gto_ordering + " format is not supported.");

}



