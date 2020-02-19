#include <iostream>
#include <vector>
#include "overlap.h"
#include "BasisSet.h"
#include <armadillo>
#include "utils.h"
#include "transforms.h"
#include "System.h"
#include "numerical.h"
using namespace std;

int main(int argc, char **argv)
{
	System my_system("N2.xyz","bas.bas");
	arma::mat Smat(my_system.bs.num_carts(),my_system.bs.num_carts());
	compute_analytical_overlap(my_system.bs,Smat);
	//uniform_cart_norm(Smat,my_system.bs);
	Smat.print();
	std::cout << std::endl;
	arma::mat spherical_ints(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(Smat,spherical_ints,my_system.bs);
	spherical_ints.print();
	arma::mat Smat2(my_system.bs.num_carts(),my_system.bs.num_carts());
	compute_numerical_overlap_mat(Smat2, my_system.bs, my_system.atoms);
	std::cout << std::endl;
	Smat2.print();
	return 0;
}

