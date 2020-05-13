/*
 * ProjectedCAP.cpp
 *
 *  Created on: May 13, 2020
 *      Author: JG
 */

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


Projected_CAP::Projected_CAP(System my_sys)
{
	//nstates = 0;
	system = my_sys;
	parameters=system.parameters;
	try
	{
		stringstream ss(parameters["nstates"]);
		ss >> nstates;
		verify_method();
		read_in_zero_order_H();
		read_in_dms();
	}
	catch(exception &e)
	{
		opencap_rethrow("Failed to initialize Projected CAP calculation.");
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
	std::cout << "Successfully read in zeroth order Hamiltonian from file:" << parameters["h0_file"] << std::endl;
	return h0;
}

void Projected_CAP::read_in_zero_order_H()
{
	if (parameters.find("h0_file")!=parameters.end())
		ZERO_ORDER_H = read_h0_file();
	else if (parameters["package"]=="qchem")
		ZERO_ORDER_H = read_qchem_energies(nstates,parameters["method"],parameters["qc_output"]);
	else if (parameters["package"]=="openmolcas")
		ZERO_ORDER_H = read_mscaspt2_heff(nstates,parameters["molcas_output"]);
	else
		std::cout << "Only q-chem and openmolcas formats are supported." << std::endl;
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
		}
		catch (exception &e)
		{
			opencap_rethrow("Failed to read TDMs from rassi.h5.");
		}
	}
	else
		opencap_throw("Only OpenMolcas rassi.h5 and Q-Chem .fchk are supported.");

}

void Projected_CAP::compute_cap_correlated_basis()
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

void Projected_CAP::reorder_cap()
{
	std::string pkg = parameters["package"];
	if (pkg=="q-chem"||pkg=="qchem")
		to_molden_ordering(AO_CAP_MAT, system.bs);
	else if(pkg=="openmolcas")
		to_molcas_ordering(AO_CAP_MAT,system.bs,system.atoms);
	else
		opencap_throw("Error. Package: " + pkg + " is unsupported.");
}

void Projected_CAP::compute_cap_matrix()
{
	check_overlap_matrix();
	std::cout << "Verified overlap matrix." << std::endl;
	std::cout << "Calculating CAP matrix in AO basis..." << std::endl;
	CAP cap_integrator(system.atoms,parameters);
	arma::mat cap_mat(system.bs.num_carts(),system.bs.num_carts());
	cap_mat.zeros();
	auto start = std::chrono::high_resolution_clock::now();
	cap_integrator.compute_cap_mat(cap_mat,system.bs);
	auto stop = std::chrono::high_resolution_clock::now();
	std::cout << "Integration time:" << std::chrono::duration<double>(stop-start).count() << std::endl;
	uniform_cart_norm(cap_mat,system.bs);
	arma::mat cap_spherical(system.bs.Nbasis,system.bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,system.bs);
	AO_CAP_MAT = cap_spherical;
	//re-order cap matrix based on electronic structure package dms came from
	reorder_cap();
	compute_cap_correlated_basis();
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
				if (abs(overlap_mat(i,j)-OVERLAP_MAT(i,j))>1E-12)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
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


void Projected_CAP::verify_method()
{
	std::string package_name = parameters["package"];
	if(parameters.find("method")==parameters.end())
		opencap_throw("Error: missing the 'method' keyword. Please choose a supported package/method.");
	std::string method = parameters["method"];
	if (package_name=="qchem")
	{
		std::vector<std::string> supported = {"eomea","eomee","eomip"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
			opencap_throw("Error: unsupported Q-Chem method. OpenCAP currently supports: 'eomea','eomee','eomip'.");
		if (parameters.find("fchk_file")==parameters.end())
			opencap_throw("Error: missing keyword: fchk_file.");
		if (parameters.find("qc_output")==parameters.end() && parameters.find("h0_file")==parameters.end())
			opencap_throw("Error: Need to specify zeroth order Hamiltonian via \"qc_output\" or \"h0_file\" fields.");
	}
	else if(package_name=="openmolcas")
	{
		std::vector<std::string> supported = {"ms-caspt2","xms-caspt2","pc-nevpt2","sc-nevpt2"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
			opencap_throw("Error: unsupported OpenMolcas method. OpenCAP currently supports: 'ms-caspt2','xms-caspt2'");
		if (parameters.find("rassi_h5")==parameters.end())
			opencap_throw("Error: missing keyword: rassi_h5.");
		if (parameters.find("molcas_output")==parameters.end() && parameters.find("h0_file")==parameters.end())
			opencap_throw("Error: Need to specify H0 Hamiltonian via \"molcas_output\" or \"h0_file\" fields.");
	}
	else
		opencap_throw("Error: unsupported package. Only QChem and OpenMolcas are currently supported.");
}
