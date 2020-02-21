#include <iostream>
#include <vector>
#include "overlap.h"
#include "BasisSet.h"
#include <armadillo>
#include "utils.h"
#include "transforms.h"
#include "System.h"
#include "numerical.h"
#include "CAP.h"
using namespace std;

int main(int argc, char **argv)
{

	System my_system("N2.xyz","bas.bas");
	std::cout << voronoi_cap(2.76,my_system.atoms).eval_pot(-0.428817,-4.094991,-2.947277) << std::endl;
	arma::mat Smat(my_system.bs.num_carts(),my_system.bs.num_carts());
	compute_analytical_overlap(my_system.bs,Smat);
	uniform_cart_norm(Smat,my_system.bs);
	arma::mat spherical_ints(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(Smat,spherical_ints,my_system.bs);
	spherical_ints.print();
	arma::mat Smat2(my_system.bs.num_carts(),my_system.bs.num_carts());
	Smat2.zeros();
	compute_numerical_overlap_mat(Smat2, my_system.bs, my_system.atoms);
	uniform_cart_norm(Smat2,my_system.bs);
	//Smat2.print();
	std::cout << std::endl;
	arma::mat spherical_ints2(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(Smat2,spherical_ints2,my_system.bs);
	spherical_ints2.print();
	/*
	arma::mat cap_mat(my_system.bs.num_carts(),my_system.bs.num_carts());
	cap_mat.zeros();
	auto cap = box_cap(2.76,2.76,4.880);
	compute_cap_mat(cap_mat, my_system.bs, my_system.atoms, cap);
    Smat2.print();
	*/

	return 0;
}

