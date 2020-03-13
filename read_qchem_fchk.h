/*
 * read_qchem_fchk.h
 *
 *  Created on: Mar 12, 2020
 *      Author: JG
 */

#ifndef READ_QCHEM_FCHK_H_
#define READ_QCHEM_FCHK_H_
#include <armadillo>

std::vector<std::string> split(const std::string& s, char delimiter);
std::array<std::vector<std::vector<arma::mat>>,2> qchem_read_in_dms(std::string dmat_filename,size_t nstates, size_t num_bf);


#endif /* READ_QCHEM_FCHK_H_ */
