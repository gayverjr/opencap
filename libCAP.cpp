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
#include "H5Cpp.h"
#include "molden_transform.h"
#include "read_qchem_fchk.h"
using namespace std;

int main(int argc, char **argv)
{
	size_t nstates = 3;
	System my_system("N2.xyz","bas.bas");
	//overlap matrix, for testing purposes
	/*
	arma::mat Smat(my_system.bs.num_carts(),my_system.bs.num_carts());
	std::cout << my_system.bs.Nbasis << std::endl;
	compute_analytical_overlap(my_system.bs,Smat);
	uniform_cart_norm(Smat,my_system.bs);
	arma::mat spherical_ints(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(Smat,spherical_ints,my_system.bs);
	to_molden_ordering(spherical_ints, my_system.bs);
	spherical_ints.print();
	*/

	//build cap matrix
	arma::mat cap_mat(my_system.bs.num_carts(),my_system.bs.num_carts());
	cap_mat.zeros();
	compute_numerical_overlap_mat(cap_mat, my_system.bs, my_system.atoms);
	uniform_cart_norm(cap_mat,my_system.bs);
	arma::mat cap_spherical(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,my_system.bs);
	to_molden_ordering(cap_spherical, my_system.bs);
	//cap_spherical.print();

	//read in DMs
	auto dms = qchem_read_in_dms("cc_test.fchk",nstates,my_system.bs.Nbasis);

	arma::mat EOMCAP(nstates,nstates);
	for (size_t row_idx=0;row_idx<EOMCAP.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<EOMCAP.n_cols;col_idx++)
		{
			EOMCAP(row_idx,col_idx) = arma::trace(dms[0][row_idx][col_idx]*cap_spherical)+
									  arma::trace(dms[1][row_idx][col_idx]*cap_spherical);
		}
	}
	std::cout << "printing eom cap" << std::endl;
	EOMCAP.print();
	/*
	arma::mat Smat2(my_system.bs.num_carts(),my_system.bs.num_carts());
	Smat2.zeros();
	compute_numerical_overlap_mat(Smat2, my_system.bs, my_system.atoms);
	uniform_cart_norm(Smat2,my_system.bs);
	//Smat2.print();
	std::cout << std::endl;
	arma::mat spherical_ints2(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(Smat2,spherical_ints2,my_system.bs);
	spherical_ints2.print();
	arma::mat cap_mat(my_system.bs.num_carts(),my_system.bs.num_carts());
	cap_mat.zeros();
	auto cap = box_cap(2.76,2.76,4.880);
	compute_cap_mat(cap_mat, my_system.bs, my_system.atoms, cap);
    Smat2.print();
	*/

	return 0;
}

