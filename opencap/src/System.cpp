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
#include <armadillo>
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
#include <carma/carma.h>

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
		arma::mat Smat(bs.num_carts(),bs.num_carts());
		compute_analytical_overlap(bs,Smat);
		uniform_cart_norm(Smat,bs);
		arma::mat spherical_ints(bs.Nbasis,bs.Nbasis);
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
	arma::mat Smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,Smat);
	uniform_cart_norm(Smat,bs);
	arma::mat spherical_ints(bs.Nbasis,bs.Nbasis);
	cart2spherical(Smat,spherical_ints,bs);
	OVERLAP_MAT = spherical_ints;
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

py::array System::get_overlap_mat()
{
	return carma::mat_to_arr(OVERLAP_MAT);
}



