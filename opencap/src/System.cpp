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
#include "utils.h"
#include "transforms.h"
#include "gto_ordering.h"
#include "CAP.h"
#include "overlap.h"
#include <cmath>
#include <limits>
#include "opencap_exception.h"

System::System(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	try
	{
		verify_system_parameters(params);
		parameters=params;
		atoms = geometry;
		if(parameters["bohr_coordinates"]!="true")
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
