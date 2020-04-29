/*
 * read_qchem_fchk.h
 *
 *  Created on: Mar 12, 2020
 *      Author: JG
 */

#ifndef QCHEM_INTERFACE_H_
#define QCHEM_INTERFACE_H_
#include <armadillo>

std::vector<std::string> split(const std::string& s, char delimiter);
std::array<std::vector<std::vector<arma::mat>>,2> qchem_read_in_dms(std::string dmat_filename,size_t nstates, size_t num_bf);
arma::mat qchem_read_overlap(std::string dmat_filename, size_t num_bf);
arma::mat read_qchem_energies(size_t nstates,std::string method,std::string output_file);

#endif /* QCHEM_INTERFACE_H_*/
