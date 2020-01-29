#include <iostream>
#include <vector>
#include <math.h>


double gauss_integral(double exp_a,double exp_b, std::vector<int> shell_a,
		std::vector<int>shell_b, std::vector<double> coord_a, std::vector<double> coord_b);
double overlap_integral(BasisFunction a, BasisFunction b);
double mcmurchie_davidson(double exp_a, double exp_b, int t, int angmom_a, int angmom_b, double dist);
