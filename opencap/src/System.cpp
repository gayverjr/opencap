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
#include "qchem_interface.h"
#include "utils.h"
#include "transforms.h"
#include "gto_ordering.h"
#include "molcas_interface.h"
#include "CAP.h"
#include "overlap.h"
#include <cmath>
#include <limits>
#include "opencap_exception.h"


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
	else if(pkg=="openmolcas")
		to_molcas_ordering(AO_CAP_MAT,bs,atoms);
	else
		opencap_throw("Error. Package: " + pkg + " is unsupported.");
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
	std::cout << "Integration time:" << std::chrono::duration<double>(stop-start).count() << std::endl;
	uniform_cart_norm(cap_mat,bs);
	arma::mat cap_spherical(bs.Nbasis,bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,bs);
	AO_CAP_MAT = cap_spherical;
	//re-order cap matrix based on electronic structure package dms came from
	reorder_cap();
	compute_cap_correlated_basis();
}

void System::check_overlap_matrix()
{
	//get overlap matrix
	arma::mat Smat(bs.num_carts(),bs.num_carts());
	compute_analytical_overlap(bs,Smat);
	uniform_cart_norm(Smat,bs);
	arma::mat spherical_ints(bs.Nbasis,bs.Nbasis);
	cart2spherical(Smat,spherical_ints,bs);
	OVERLAP_MAT = spherical_ints;
	if (parameters["package"]=="qchem")
	{
		to_molden_ordering(spherical_ints, bs);
		auto qchem_smat = qchem_read_overlap(parameters["fchk_file"],bs.Nbasis);
		if(OVERLAP_MAT.n_rows != qchem_smat.n_rows || OVERLAP_MAT.n_cols != qchem_smat.n_cols)
			opencap_throw("Basis set has wrong dimension when checking overlap matrix. Check your input files.");
		bool conflicts = false;
		for (size_t i=0;i<qchem_smat.n_rows;i++)
		{
			for(size_t j=0;j<qchem_smat.n_cols;j++)
			{
				if (abs(qchem_smat(i,j)-spherical_ints(i,j))>1E-5)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
					std::cout << "Qchem says:" << qchem_smat(i,j) << std::endl;
					std::cout << "OpenCAP says:" << spherical_ints(i,j) << std::endl;
					conflicts = true;
				}
			}
		}
			if(conflicts)
				opencap_throw("Error. Overlap matrices do not match. Check your basis set.");
	}
	else if (parameters["package"]=="openmolcas")
	{
		to_molcas_ordering(spherical_ints,bs,atoms);
		arma::mat overlap_mat = read_rassi_overlap(parameters["rassi_h5"]);
		if(OVERLAP_MAT.n_rows != overlap_mat.n_rows || OVERLAP_MAT.n_cols != overlap_mat.n_cols)
			opencap_throw("Basis set has wrong dimension when checking overlap matrix. Check your input files.");
		std::cout << std::fixed << std::setprecision(10);
		bool conflicts = false;
		for (size_t i=0;i<overlap_mat.n_rows;i++)
		{
			for(size_t j=0;j<overlap_mat.n_cols;j++)
			{
				if (abs(overlap_mat(i,j)-spherical_ints(i,j))>1E-12)
				{
					std::cout << "Conflict at:" << i << "," << j << std::endl;
					std::cout << "OpenMolcas says:" << overlap_mat(i,j) << std::endl;
					std::cout << "OpenCAP says:" << spherical_ints(i,j) << std::endl;
					conflicts = true;
				}
			}
		}
		if (conflicts)
			opencap_throw("Error. Overlap matrices do not match. Check your basis set.");
	}
	else
		opencap_throw("Error. Package: " + parameters["package"] + " is unsupported.");
}

void System::read_in_dms()
{

	stringstream ss(parameters["nstates"]);
	size_t nstates;
	ss >> nstates;
	if (parameters["package"]=="qchem")
	{
		std::string dmat_filename = parameters["fchk_file"];
		try
		{
			auto parsed_dms = qchem_read_in_dms(dmat_filename,nstates,bs.Nbasis);
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

arma::mat System::read_h0_file()
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
	return h0;
}

void System::read_in_zero_order_H()
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

System::System(std::vector<Atom> geometry,std::map<std::string, std::string> params)
{
	parameters=params;
	atoms = geometry;
	if(parameters["bohr_coordinates"]=="false")
	{
		for (size_t i=0;i<atoms.size();i++)
			atoms[i].ang_to_bohr();
	}
	bs = BasisSet(atoms,parameters);
	try
	{
		check_overlap_matrix();
		std::cout << "Verified overlap matrix." << std::endl;
		stringstream ss(parameters["nstates"]);
		ss >> nstates;
		read_in_dms();
		read_in_zero_order_H();
		std::cout << "Successfully read in zeroth order Hamiltonian." << std::endl;
	}
	catch (exception& e)
	{
		opencap_rethrow("Failed to construct system.");
	}

}
