/*
 * utils.h
 *
 *  Created on: Feb 12, 2020
 *      Author: JG
 */
#include "Shell.h"
#include <armadillo>
double fact2(int n);
double fact(int n);
std::vector<std::array<size_t,3>> libcap_carts_ordering(Shell shell);
std::vector<int> libcap_harmonic_ordering(Shell shell);
int parity(int i);
double binom(int p, int q);
std::vector<std::string> split(const std::string& s, char delimiter);
void fill_LT(std::vector<double> matrix_elements, arma::mat &opdm);
void fill_mat(std::vector<double> &matrix_elements, arma::mat &opdm);
