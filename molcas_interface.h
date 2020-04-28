/*
 * readMolcasHDF5.h
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */
#include <armadillo>
#ifndef MOLCAS_INTERFACE_H_
#define MOLCAS_INTERFACE_H_

std::array<std::vector<std::vector<arma::mat>>,2> read_rassi_tdms(std::string dmat_filename);
arma::mat read_rassi_overlap(std::string filename);
arma::mat read_mscaspt2_heff(size_t nstates, std::string filename);

#endif /* MOLCAS_INTERFACE_H_ */
