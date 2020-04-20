/*
 * System.h
 *
 *  Created on: Feb 18, 2020
 *      Author: JG
 */
#include "Atom.h"
#include "BasisSet.h"
#include "CAP.h"
#include <map>
#include <armadillo>
#ifndef SYSTEM_H_
#define SYSTEM_H_

class System
{
public:
	std::vector<Atom> atoms;
	std::map<std::string, std::string> parameters;
	//std::string cart;
	BasisSet bs;
	//int nelectron;
	//System(std::string xyz_name,std::string basis_name);
	System();
	System(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	arma::mat AO_CAP_MAT;
	arma::mat CORRELATED_CAP_MAT;
	arma::mat ZERO_ORDER_H;
	std::vector<std::vector<arma::mat>> alpha_dms;
	std::vector<std::vector<arma::mat>> beta_dms;
	size_t nstates;
	void compute_cap_matrix();

private:
	bool check_overlap_matrix();
	void ang_to_bohr();
	bool verify_cap_parameters(std::string key);
	bool check_parameters();
	bool read_in_dms();
	void compute_cap_correlated_basis();
	void reorder_cap();
	bool verify_method(std::string key);
	bool read_in_zero_order_H();
	bool read_qchem_energies();
	bool read_molcas_Heff();

};



#endif /* SYSTEM_H_ */
