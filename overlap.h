#include <iostream>
#include <vector>
#include <math.h>
#include "Shell.h"
#include <armadillo>


double gauss_integral(double exp_a,std::array<size_t, 3> cart_a,std::array<double,3> coord_a,
		double exp_b,std::array<size_t, 3> cart_b, std::array<double,3> coord_b);
double overlap_integral(Shell a, std::array<size_t,3> cart_a, Shell b, std::array<size_t,3> cart_);
double mcmurchie_davidson(double exp_a, double exp_b, int t, int angmom_a, int angmom_b, double dist);
void shell_overlap(Shell shell_a, Shell shell_b,arma::subview<double>&sub_mat);
std::vector<std::array<size_t,3>> get_carts_ordering(Shell shell);
