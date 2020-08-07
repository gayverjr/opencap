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
#include <eigen3/Eigen/Dense>

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
		if(check_keyword(key,"system")||key=="geometry")
			parameters[key]=value;
		else
			opencap_throw("Invalid key:" +key);
    }
    if (parameters.find("molecule")==parameters.end()||parameters.find("basis_file")==parameters.end())
    	opencap_throw("Error: molecule and basis_file keywords are required.");
	if(compare_strings(parameters["molecule"],"qchem_fchk"))
	{
		atoms = read_geometry_from_fchk(parameters["basis_file"]);
		bs = read_basis_from_fchk(parameters["basis_file"],atoms);
	}
	else if (compare_strings(parameters["molecule"],"molcas_rassi"))
	{
		atoms = read_geometry_from_rassi(parameters["basis_file"]);
		bs = read_basis_from_rassi(parameters["basis_file"],atoms);
	}
	else if(compare_strings(parameters["molecule"],"molden"))
	{
		atoms = read_geometry_from_molden(parameters["basis_file"]);
		bs = read_basis_from_molden(parameters["basis_file"],atoms);
	}
	else if(compare_strings(parameters["molecule"],"inline"))
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
	parameters["molecule"]   = file_type;
	parameters["basis_file"] = filename;
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

void System::renormalize_overlap(Eigen::MatrixXd smat)
{
	Eigen::MatrixXd overlap_copy = OVERLAP_MAT;
	std::vector<double> scalars;
	for (size_t i=0;i<smat.rows();i++)
		scalars.push_back(sqrt(smat(i,i)));
	for (size_t i=0;i<overlap_copy.rows();i++)
	{
		for(size_t j=0;j<overlap_copy.cols();j++)
			overlap_copy(i,j)= OVERLAP_MAT(i,j)* scalars[i] * scalars[j];
	}
	for (size_t i=0;i<smat.rows();i++)
	{
		for(size_t j=0;j<smat.cols();j++)
		{
			if (abs(smat(i,j)-overlap_copy(i,j))>1E-5)
			{
				opencap_throw("Error: Could not verify overlap matrix after re-normalization."
								"Verify that your basis is specified properly, or use a different type of input. If the "
								"issue persists, please file a bug report.");
			}
		}
	}

}

bool System::check_overlap_mat(Eigen::MatrixXd smat, std::string ordering, std::string basis_file)
{
	std::vector<bf_id> ids;
	if(compare_strings(ordering,"pyscf"))
		ids = get_pyscf_ids(bs);
	else if(compare_strings(ordering,"openmolcas"))
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(bs,basis_file);
	}
	else if(compare_strings(ordering,"qchem"))
		ids = get_qchem_ids(bs);
	else if(compare_strings(ordering,"opencap"))
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
				conflicts = true;
				if( (abs(smat(i,j))<1E-10 && abs(OVERLAP_MAT(i,j))>1E-10) || (abs(smat(i,j))>1E-10 && abs(OVERLAP_MAT(i,j))<1E-10) )
				{
					opencap_throw("Error: The dimensions of the overlap matrices match, but the elements do not. "
							"Verify that your basis is specified properly, or use a different type of input. If the "
							"issue persists, please file a bug report.");
				}

			}
		}
	}
	std::string message;
	if (conflicts)
	{
		if(python)
			py::print("Warning: the overlap matrices differ numerically, but there are no non-matching "
				 "zeroes.\nIf you are using cartesian GTOs, this is expected.");
		else
			std::cout << "Warning: the overlap matrices differ numerically, but there are no non-matching "
			 "zeroes.\nIf you are using cartesian GTOs, this is expected." << std::endl;
	}
	else
	{
		message= "Verified overlap matrix.";
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
		return true;
	}
	renormalize_overlap(smat);
	if(python)
		py::print("Verified overlap matrix after re-normalization.\nPlease use the renormalize or"
				"renormalize_cap functions before computing the projected CAP.");
	else
		std::cout << "Verified overlap matrix after re-normalization." << std::endl;
	return false;

}


