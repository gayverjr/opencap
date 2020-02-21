/*
 * numerical.h
 *
 *  Created on: Feb 18, 2020
 *      Author: JG
 */
#include "BasisSet.h"
#include <armadillo>



void compute_numerical_overlap_mat(arma::mat &Smat, BasisSet bs,std::vector<Atom> atoms);
void compute_cap_mat(arma::mat &Smat, BasisSet bs,std::vector<Atom> atoms, cap_pot &pot);


