#include <functional>
#include <Eigen/Dense>
#include "BasisSet.h"
#include "Atom.h"


void integrate_cap_numerical(Eigen::MatrixXd &cap_mat, BasisSet bs, std::vector<Atom> atoms, 
double radial_precision, size_t angular_points,
std::function<std::vector<double>(std::vector<double> &, std::vector<double> &, 
		std::vector<double> &, std::vector<double> &, int)> &cap_func);