/*
 * readMolcasHDF5.h
 *
 *  Created on: Mar 15, 2020
 *      Author: JG
 */
#include <armadillo>
#ifndef READ_RASSI_h5_H_
#define READ_RASSI_h5_H_

std::array<std::vector<std::vector<arma::mat>>,2> read_rassi_tdms(std::string dmat_filename);
arma::mat read_rassi_overlap(std::string filename);
arma::mat read_mscaspt2_heff(size_t nstates, std::string filename);

#endif /* READMOLCASHDF5_H_ */
