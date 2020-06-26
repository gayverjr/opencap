/*
 * ProjectedCAP.cpp
 *
 *  Created on: May 13, 2020
 *      Author: JG
 */
#include <pybind11/pybind11.h>
#include "ProjectedCAP.h"
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
#include "qchem_interface.h"
#include "molcas_interface.h"
#include "CAP.h"
#include "overlap.h"
#include <cmath>
#include <limits>
#include "opencap_exception.h"
#include <carma/carma.h>

Projected_CAP::Projected_CAP(System my_sys,std::map<std::string, std::string> params)
{
	system = my_sys;
	try
	{
		verify_method(params);
		parameters=params;
		stringstream ss(parameters["nstates"]);
		ss >> nstates;
		read_in_zero_order_H();
		read_in_dms();
	}
	catch(exception &e)
	{
		opencap_rethrow("Failed to initialize Projected CAP calculation.");
	}
	python = false;
}

Projected_CAP::Projected_CAP(System my_sys, py::dict dict, size_t num_states, std::string gto_ordering)
{
	std::vector<std::string> valid_keywords = {"cap_type","cap_x","cap_y","cap_z",
			"r_cut","radial_precision","angular_points"};
	std::map<std::string, std::string> params;
    for (auto item : dict)
    {
    	std::string key = py::str(item.first).cast<std::string>();
    	std::string value = py::str(item.second).cast<std::string>();
		transform(key.begin(),key.end(),key.begin(),::tolower);
		transform(value.begin(),value.end(),value.begin(),::tolower);
    	if (std::find(valid_keywords.begin(),
    			valid_keywords.end(),key)==valid_keywords.end())
    		opencap_throw("Invalid key in dictionary:`" + key + "'\n");
		params[key]=value;
    }
    try
    {
    	CAP cap_integrator(system.atoms,params);
		for (auto item: params)
			parameters[item.first]=item.second;
		system = my_sys;
		nstates = num_states;
		if(num_states<1)
			opencap_throw("Error: not enough states to run calculation.");
		std::transform(gto_ordering.begin(), gto_ordering.end(), gto_ordering.begin(), ::tolower);
		if(!(gto_ordering=="openmolcas"||gto_ordering=="qchem" || gto_ordering=="pyscf"))
			opencap_throw("Error: " + gto_ordering + " ordering is unsupported.");
		parameters["package"] = gto_ordering;
		alpha_dms = std::vector<std::vector<arma::mat>>(nstates,
				std::vector<arma::mat>(nstates));
		beta_dms = std::vector<std::vector<arma::mat>>(nstates,
				std::vector<arma::mat>(nstates));
		python = true;
    }
    catch(exception &e)
    {
    	opencap_throw("Failed to construct Projected CAP object.");
    }

}

arma::mat Projected_CAP::read_h0_file()
{
	arma::mat h0(nstates,nstates);
	h0.zeros();
	std::ifstream is(parameters["h0_file"]);
	if (is.good())
	{
		//first line should be diagonal or full
		std::string line;
		std::getline(is,line);
		std::string mat_type = line;
		std::transform(mat_type.begin(), mat_type.end(), mat_type.begin(), ::tolower);
		if (mat_type=="diagonal")
		{
			for(size_t i=0;i<nstates;i++)
			{
				std::getline(is,line);
				h0(i,i)=std::stod(line);
			}
		}
		else if (mat_type=="full")
		{
			for (size_t i=0;i<nstates;i++)
			{
				std::getline(is,line);
				std::vector<std::string> tokens = split(line,' ');
				for(size_t j=0;j<nstates;j++)
				{
					h0(i,j)=std::stod(tokens[j]);
				}
			}
		}
		else
		{
			return h0;
		}
	}
	std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["h0_file"];
	if(python)
		py::print(message);
	else
		std::cout << message << std::endl;
	return h0;
}

void Projected_CAP::read_in_zero_order_H()
{
	if (parameters.find("h0_file")!=parameters.end())
		ZERO_ORDER_H = read_h0_file();
	else if (parameters["package"]=="qchem")
	{
		ZERO_ORDER_H = read_qchem_energies(nstates,parameters["method"],parameters["qc_output"]);
		std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["qc_output"];
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
	}
	else if (parameters["package"]=="openmolcas")
	{
		ZERO_ORDER_H = read_mscaspt2_heff(nstates,parameters["molcas_output"]);
		std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["molcas_output"];
		if(python)
			py::print(message);
		else
			std::cout << message << std::endl;
	}
	else
		opencap_throw("Error: Only q-chem and openmolcas formats are supported.");
}

void Projected_CAP::read_in_dms()
{
	if (parameters["package"]=="qchem")
	{
		std::string dmat_filename = parameters["fchk_file"];
		try
		{
			auto parsed_dms = qchem_read_in_dms(dmat_filename,nstates,system.bs.Nbasis);
			alpha_dms = parsed_dms[0];
			beta_dms = parsed_dms[1];
			std::string message= "Successfully read in densities from file:" + parameters["fchk_file"];
			if(python)
				py::print(message);
			else
				std::cout << message << std::endl;
		}
		catch (exception &e)
		{
			opencap_rethrow("Failed to read TDMs from Q-Chem .fchk.");
		}
	}
	else if (parameters["package"]=="openmolcas")
	{
		try
		{
			auto parsed_dms = read_rassi_tdms(parameters["rassi_h5"]);
			alpha_dms = parsed_dms[0];
			beta_dms = parsed_dms[1];
			std::string message = "Successfully read in densities from file:" + parameters["rassi_h5"];
			if(python)
				py::print(message);
			else
				std::cout << message << std::endl;
		}
		catch (exception &e)
		{
			opencap_rethrow("Failed to read TDMs from file:"+parameters["rassi_h5"]);
		}
	}
	else
		opencap_throw("Only OpenMolcas rassi.h5 and Q-Chem .fchk are supported.");

}

void Projected_CAP::compute_projected_cap()
{
	verify_data();
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

void Projected_CAP::reorder_cap()
{
	std::string pkg = parameters["package"];
	if (pkg=="q-chem"||pkg=="qchem")
		to_molden_ordering(AO_CAP_MAT, system.bs);
	else if(pkg=="openmolcas")
		to_molcas_ordering(AO_CAP_MAT,system.bs,system.atoms);
	else if(pkg == "pyscf")
		to_pyscf_ordering(AO_CAP_MAT,system.bs);
	else
		opencap_throw("Error. Package: " + pkg + " is unsupported.");
}

void Projected_CAP::compute_ao_cap()
{
	if(!python)
	{
		check_overlap_matrix();
		std::cout << "Verified overlap matrix." << std::endl;
	}
	if(python)
		py::print("Calculating CAP matrix in AO basis...");
	else
		std::cout << "Calculating CAP matrix in AO basis..." << std::endl;
	CAP cap_integrator(system.atoms,parameters);
	arma::mat cap_mat(system.bs.num_carts(),system.bs.num_carts());
	cap_mat.zeros();
	auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_cap_mat(cap_mat,system.bs);
	auto stop = std::chrono::high_resolution_clock::now();
	auto total_time = std::chrono::duration<double>(stop-start).count();
	if(python)
		py::print("Integration time:"+std::to_string(total_time));
	else
		std::cout << "Integration time:" << std::to_string(total_time) << std::endl;
	uniform_cart_norm(cap_mat,system.bs);
	arma::mat cap_spherical(system.bs.Nbasis,system.bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,system.bs);
	AO_CAP_MAT = cap_spherical;
	//re-order cap matrix based on electronic structure package dms came from
	reorder_cap();
}

void Projected_CAP::check_overlap_matrix()
{
	//get overlap matrix
	arma::mat OVERLAP_MAT=system.OVERLAP_MAT;
	if (parameters["package"]=="qchem")
	{
		to_molden_ordering(OVERLAP_MAT, system.bs);
		auto qchem_smat = qchem_read_overlap(parameters["fchk_file"],system.bs.Nbasis);
		if(OVERLAP_MAT.n_rows != qchem_smat.n_rows || OVERLAP_MAT.n_cols != qchem_smat.n_cols)
			opencap_throw("Basis set has wrong dimension when checking overlap matrix. "
					"Verify that cart_bf field is set correctly (default is all spherical harmonic), "
					"and that the basis set is the same as that used in the electronic structure"
					" calculation.");
		bool conflicts = false;
		for (size_t i=0;i<qchem_smat.n_rows;i++)
		{
			for(size_t j=0;j<qchem_smat.n_cols;j++)
			{
				if (abs(qchem_smat(i,j)-OVERLAP_MAT(i,j))>1E-5)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
					std::cout << "Qchem says:" << qchem_smat(i,j) << std::endl;
					std::cout << "OpenCAP says:" << OVERLAP_MAT(i,j) << std::endl;
					conflicts = true;
				}
			}
		}
		if (conflicts)
			opencap_throw("Error. The dimensions of the overlap matrices match, but the elements do not. Verify that"
					" your geometry has the right units, and that your basis file is ordered the same as in "
					"the electronic structure calculation.");
	}
	else if (parameters["package"]=="openmolcas")
	{
		to_molcas_ordering(OVERLAP_MAT,system.bs,system.atoms);
		arma::mat overlap_mat = read_rassi_overlap(parameters["rassi_h5"]);
		if(OVERLAP_MAT.n_rows != overlap_mat.n_rows || OVERLAP_MAT.n_cols != overlap_mat.n_cols)
			opencap_throw("Basis set has wrong dimension when checking overlap matrix. "
					"Verify that cart_bf field is set correctly (default is all spherical harmonic), "
					"and that the basis set is the same as that used in the electronic structure"
					" calculation.");
		std::cout << std::fixed << std::setprecision(10);
		bool conflicts = false;
		for (size_t i=0;i<overlap_mat.n_rows;i++)
		{
			for(size_t j=0;j<overlap_mat.n_cols;j++)
			{
				if (abs(overlap_mat(i,j)-OVERLAP_MAT(i,j))>1E-5)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
					std::cout <<  overlap_mat(i,j)-OVERLAP_MAT(i,j) << std::endl;
					std::cout << "OpenMolcas says:" << overlap_mat(i,j) << std::endl;
					std::cout << "OpenCAP says:" << OVERLAP_MAT(i,j) << std::endl;
					conflicts = true;
				}
			}
		}
		if (conflicts)
			opencap_throw("Error. The dimensions of the overlap matrices match, but the elements do not. Verify that"
					" your geometry has the right units, and that your basis file is ordered the same as in "
					"the electronic structure calculation.");
	}
	else
		opencap_throw("Error. Package: " + parameters["package"] + " is unsupported.");
}


void Projected_CAP::verify_method(std::map<std::string,std::string> params)
{
	std::string package_name;
	if (!python)
	{
		if(params.find("nstates")==params.end())
			opencap_throw("Error: missing the 'nstates' keyword. "
					"Please specify the number of electronic states.")
		package_name = params["package"];
	}
	else
		package_name = parameters["package"];
	if(params.find("method")==params.end())
		opencap_throw("Error: missing the 'method' keyword. "
				"Please choose a supported package/method.");
	std::string method = params["method"];
	if (package_name=="qchem")
	{
		std::vector<std::string> supported = {"eomea","eomee","eomip"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
			opencap_throw("Error: unsupported Q-Chem method. OpenCAP currently supports: 'eomea','eomee','eomip'.");
		if (params.find("fchk_file")==params.end())
			opencap_throw("Error: missing keyword: fchk_file.");
		if (params.find("qc_output")==params.end() && params.find("h0_file")==params.end())
			opencap_throw("Error: Need to specify zeroth order Hamiltonian via \"qc_output\" or \"h0_file\" fields.");
	}
	else if(package_name=="openmolcas")
	{
		std::vector<std::string> supported = {"ms-caspt2","xms-caspt2","pc-nevpt2","sc-nevpt2"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
			opencap_throw("Error: unsupported OpenMolcas method. OpenCAP currently supports: 'ms-caspt2','xms-caspt2'");
		if (params.find("rassi_h5")==params.end())
			opencap_throw("Error: missing keyword: rassi_h5.");
		if (params.find("molcas_output")==params.end() && params.find("h0_file")==params.end())
			opencap_throw("Error: Need to specify H0 Hamiltonian via \"molcas_output\" or \"h0_file\" fields.");
	}
	else
		opencap_throw("Error: unsupported package. Only QChem and OpenMolcas are currently supported.");
}

void Projected_CAP::run()
{
	try
	{
		compute_ao_cap();
		compute_projected_cap();
	}
	catch(exception &e)
	{
		opencap_rethrow("Error initializing projected CAP calculation.");
	}
}

void Projected_CAP::verify_data()
{
	//check that number of states is correct
	if(!(nstates == alpha_dms.size() && nstates == beta_dms.size()))
		opencap_throw("Error: dimensions of H and density matrix container do not match number of states.")
	//check that dimensionality of density matrices is correct
	//alpha
	for(size_t i=0;i<nstates;i++)
	{
		for(size_t j=0;j<nstates;j++)
			if(!(alpha_dms[i][j].n_cols == alpha_dms[i][j].n_rows && alpha_dms[i][j].n_rows == system.bs.Nbasis))
				opencap_throw("Error: Dimensionality of the density matrices do not match"
						"the basis set specified in the system object.")
	}
	//beta
	for(size_t i=0;i<nstates;i++)
	{
		for(size_t j=0;j<nstates;j++)
		{
			if(!(beta_dms[i][j].n_cols == beta_dms[i][j].n_rows && beta_dms[i][j].n_rows == system.bs.Nbasis))
				opencap_throw("Error: Dimensionality of the density matrices do not match"
			"the basis set specified in the system object.")
		}
	}
}

py::array Projected_CAP::get_ao_cap()
{
	return carma::mat_to_arr(AO_CAP_MAT);
}

py::array Projected_CAP::get_projected_cap()
{
	return carma::mat_to_arr(CORRELATED_CAP_MAT);
}

py::array Projected_CAP::get_H()
{
	return carma::mat_to_arr(ZERO_ORDER_H);
}

void Projected_CAP::add_tdm(py::array_t<double> & alpha_density,
		py::array_t<double> & beta_density,size_t row_idx, size_t col_idx)
{
	arma::mat alpha_dm = carma::arr_to_mat<double>(alpha_density);
	arma::mat beta_dm = carma::arr_to_mat<double>(beta_density);
	if(!(alpha_dm.n_cols == alpha_dm.n_rows && alpha_dm.n_cols == system.bs.Nbasis &&
			beta_dm.n_cols == beta_dm.n_rows && beta_dm.n_cols == system.bs.Nbasis))
		opencap_throw("Error: Dimensionality of the density matrices "
				"do not match the basis set specified in the system object.");
	alpha_dms [row_idx][col_idx] = alpha_dm;
	beta_dms [row_idx][col_idx] = beta_dm;
}

void Projected_CAP::add_tdm(py::array_t<double> & tdm,size_t row_idx, size_t col_idx)
{
	arma::mat dmat = 0.5 * carma::arr_to_mat<double>(tdm);
	if(!(dmat.n_cols == dmat.n_rows&& dmat.n_cols == system.bs.Nbasis))
		opencap_throw("Error: Dimensionality of the density matrix "
				"does not match the basis set specified in the system object.");
	alpha_dms[row_idx][col_idx] = dmat;
	beta_dms[row_idx][col_idx] = dmat;
}

void Projected_CAP::set_h0(py::array_t<double> &h0)
{
	arma::mat mat = carma::arr_to_mat<double>(h0);
	if(!(mat.n_cols==mat.n_rows&&mat.n_cols==nstates))
		opencap_throw("Error: dimensionality of matrix does not match number of states.");
	ZERO_ORDER_H = mat ;
}

void Projected_CAP::read_electronic_structure_data(py::dict dict)
{
	std::vector<std::string> valid_keywords = {"method","qc_output","h0_file","rassi_h5",
			"fchk_file","molcas_output"};
	std::map<std::string, std::string> params;
    for (auto item : dict)
    {
    	std::string key = py::str(item.first).cast<std::string>();
    	std::string value = py::str(item.second).cast<std::string>();
		transform(key.begin(),key.end(),key.begin(),::tolower);
		transform(value.begin(),value.end(),value.begin(),::tolower);
    	if (std::find(valid_keywords.begin(),
    			valid_keywords.end(),key)==valid_keywords.end())
    		opencap_throw("Invalid key in dictionary:`" + key + "'\n");
		params[key]=value;
    }
	try
	{
		verify_method(params);
		for (auto item: params)
			parameters[item.first]=item.second;
		read_in_zero_order_H();
		read_in_dms();
	}
	catch(exception &e)
	{
		//remove bad params
		for (auto item: dict)
		{
			std::string key = py::str(item.first).cast<std::string>();
	    	std::string value = py::str(item.second).cast<std::string>();
			if(parameters[key] == value)
			{
				auto it1 = parameters.find(key);
				parameters.erase(it1);
			}
		}
		opencap_rethrow("Failed to read electronic structure data.");
	}
}
