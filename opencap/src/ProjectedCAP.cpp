/*
 * ProjectedCAP.cpp
 *
 *  Created on: May 13, 2020
 *      Author: JG
 */
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "ProjectedCAP.h"
#include "System.h"
#include "BasisSet.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include "CAP.h"
#include <string>
#include <map>
#include <list>
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
#include <Eigen/Dense>


Projected_CAP::Projected_CAP(System my_sys,std::map<std::string, std::string> params)
{
	system = my_sys;
	try
	{
		python = false;
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
}

Projected_CAP::Projected_CAP(System my_sys, py::dict dict, size_t num_states, std::string gto_ordering)
{
	std::vector<std::string> valid_keywords = {"cap_type","cap_x","cap_y","cap_z",
			"r_cut","radial_precision","angular_points"};
	std::map<std::string, std::string> params;
	python = true;
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
		alpha_dms = std::vector<std::vector<Eigen::MatrixXd>>(nstates,
				std::vector<Eigen::MatrixXd>(nstates));
		beta_dms = std::vector<std::vector<Eigen::MatrixXd>>(nstates,
				std::vector<Eigen::MatrixXd>(nstates));
    }
    catch(exception &e)
    {
    	opencap_throw("Failed to construct Projected CAP object.");
    }

}

Eigen::MatrixXd Projected_CAP::read_h0_file()
{
	Eigen::MatrixXd h0(nstates,nstates);
	h0= Eigen::MatrixXd::Zero(nstates,nstates);
	std::cout << "Reading in Hamiltonian from:" << parameters["h0_file"] << std::endl;
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
		ZERO_ORDER_H = read_qchem_energies(nstates,parameters["method"],parameters["qchem_output"]);
		std::string message = "Successfully read in zeroth order Hamiltonian from file:" + parameters["qchem_output"];
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
		std::string dmat_filename = parameters["qchem_fchk"];
		try
		{
			if(parameters["method"]=="eomee")
			{
				auto parsed_dms = qchem_read_in_dms_closed_shell(dmat_filename,nstates,system.bs);
				alpha_dms = parsed_dms[0];
				beta_dms = parsed_dms[1];
			}
			else
			{
				auto parsed_dms = qchem_read_in_dms_open_shell(dmat_filename,nstates,system.bs);
				alpha_dms = parsed_dms[0];
				beta_dms = parsed_dms[1];
			}
			std::string message= "Successfully read in densities from file:" + parameters["qchem_fchk"];
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
			auto parsed_dms = read_rassi_tdms(parameters["rassi_h5"],system.bs);
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
	Eigen::MatrixXd EOMCAP(nstates,nstates);
	EOMCAP= Eigen::MatrixXd::Zero(nstates,nstates);
    std::cout << std::fixed << std::setprecision(10);
	for (size_t row_idx=0;row_idx<EOMCAP.rows();row_idx++)
	{
		for (size_t col_idx=0;col_idx<EOMCAP.cols();col_idx++)
		{
			Eigen::MatrixXd alpha_mat_prod = alpha_dms[row_idx][col_idx]*AO_CAP_MAT;
			Eigen::MatrixXd beta_mat_prod = beta_dms[row_idx][col_idx]*AO_CAP_MAT;
			EOMCAP(row_idx,col_idx) =  alpha_mat_prod.trace()+beta_mat_prod.trace();
			EOMCAP(row_idx,col_idx) = -1.0* EOMCAP(row_idx,col_idx);
		}
	}
	CORRELATED_CAP_MAT = EOMCAP;
}

void Projected_CAP::compute_ao_cap()
{
	if(python)
		py::print("Calculating CAP matrix in AO basis...");
	else
		std::cout << "Calculating CAP matrix in AO basis..." << std::endl;
	CAP cap_integrator(system.atoms,parameters);
	Eigen::MatrixXd cap_mat(system.bs.num_carts(),system.bs.num_carts());
	cap_mat= Eigen::MatrixXd::Zero(system.bs.num_carts(),system.bs.num_carts());
	auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_cap_mat(cap_mat,system.bs);
	auto stop = std::chrono::high_resolution_clock::now();
	auto total_time = std::chrono::duration<double>(stop-start).count();
	if(python)
		py::print("Integration time:"+std::to_string(total_time));
	else
		std::cout << "Integration time:" << std::to_string(total_time) << std::endl;
	uniform_cart_norm(cap_mat,system.bs);
	Eigen::MatrixXd cap_spherical(system.bs.Nbasis,system.bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,system.bs);
	AO_CAP_MAT = cap_spherical;
}

void Projected_CAP::check_overlap_matrix()
{
	//get overlap matrix
	Eigen::MatrixXd OVERLAP_MAT=system.OVERLAP_MAT;
	Eigen::MatrixXd smat;
	if (parameters["package"]=="qchem")
		smat = qchem_read_overlap(parameters["qchem_fchk"],system.bs);
	else if (parameters["package"]=="openmolcas")
		smat = read_rassi_overlap(parameters["rassi_h5"],system.bs);
	system.check_overlap_mat(smat,"opencap");
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
		if (params.find("qchem_fchk")==params.end())
			opencap_throw("Error: missing keyword: qchem_fchk.");
		if (params.find("qchem_output")==params.end() && params.find("h0_file")==params.end())
			opencap_throw("Error: Need to specify zeroth order Hamiltonian via \"qchem_output\" or \"h0_file\" fields.");
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
		check_overlap_matrix();
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
			if(!(alpha_dms[i][j].cols() == alpha_dms[i][j].rows() && alpha_dms[i][j].rows() == system.bs.Nbasis))
				opencap_throw("Error: Dimensionality of the density matrices do not match"
						"the basis set specified in the system object.")
	}
	//beta
	for(size_t i=0;i<nstates;i++)
	{
		for(size_t j=0;j<nstates;j++)
		{
			if(!(beta_dms[i][j].cols() == beta_dms[i][j].rows() && beta_dms[i][j].rows() == system.bs.Nbasis))
				opencap_throw("Error: Dimensionality of the density matrices do not match"
			"the basis set specified in the system object.")
		}
	}
}

Eigen::MatrixXd Projected_CAP::get_ao_cap()
{
	return AO_CAP_MAT;
}

Eigen::MatrixXd Projected_CAP::get_projected_cap()
{
	return CORRELATED_CAP_MAT;
}

Eigen::MatrixXd Projected_CAP::get_H()
{
	return ZERO_ORDER_H;
}

void Projected_CAP::add_tdms(Eigen::MatrixXd & alpha_density,
		Eigen::MatrixXd & beta_density,size_t row_idx, size_t col_idx,
		std::string ordering,std::string basis_file)
{
	Eigen::MatrixXd alpha_dm = alpha_density;
	Eigen::MatrixXd beta_dm = beta_density;
	if(!(alpha_dm.cols() == alpha_dm.rows() && alpha_dm.cols() == system.bs.Nbasis &&
			beta_dm.cols() == beta_dm.rows() && beta_dm.cols() == system.bs.Nbasis))
		opencap_throw("Error: Dimensionality of the density matrices "
				"do not match the basis set specified in the system object.");
	std::vector<bf_id> ids;
	if(ordering=="pyscf")
		ids = get_pyscf_ids(system.bs);
	else if(ordering=="openmolcas")
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(system.bs,basis_file);
	}
	else if(ordering=="qchem"||ordering=="molden")
		ids = get_molden_ids(system.bs);
	else
		opencap_throw("Error: " + ordering +" is unsupported.");
	to_opencap_ordering(alpha_dm,system.bs,ids);
	to_opencap_ordering(beta_dm,system.bs,ids);
	alpha_dms [row_idx][col_idx] = alpha_dm;
	beta_dms [row_idx][col_idx] = beta_dm;
}

void Projected_CAP::add_tdm(Eigen::MatrixXd tdm,size_t row_idx, size_t col_idx,std::string ordering,
		std::string basis_file)
{
	Eigen::MatrixXd dmat;
	dmat = 0.5 * tdm;
	if(!(dmat.cols() == dmat.rows()&& dmat.cols() == system.bs.Nbasis))
		opencap_throw("Error: Dimensionality of the density matrix "
				"does not match the basis set specified in the system object.");
	std::vector<bf_id> ids;
	if(ordering=="pyscf")
		ids = get_pyscf_ids(system.bs);
	else if(ordering=="openmolcas")
	{
		if(basis_file=="")
			opencap_throw("Error: OpenMolcas ordering requires a valid HDF5 file "
					"specified with the basis_file optional argument.");
		ids = get_molcas_ids(system.bs,basis_file);
	}
	else if(ordering=="qchem"||ordering=="molden")
		ids = get_molden_ids(system.bs);
	else
		opencap_throw("Error: " + ordering +" is unsupported.");
	to_opencap_ordering(dmat,system.bs,ids);
	alpha_dms[row_idx][col_idx] = dmat;
	beta_dms[row_idx][col_idx] = dmat;
}

void Projected_CAP::read_electronic_structure_data(py::dict dict)
{
	std::vector<std::string> valid_keywords = {"method","qchem_output","rassi_h5","qchem_fchk","molcas_output"};
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
