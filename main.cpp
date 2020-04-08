#include <iostream>
#include <vector>
#include "System.h"
#include <chrono>
#include <map>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>
#include "InputParser.h"
#include "Atom.h"

int main(int argc, char **argv)
{
	std::cout << "Welcome to OpenCAP!" << std::endl;
	auto t_start = std::chrono::high_resolution_clock::now();
	std::tuple<std::vector<Atom>,std::map<std::string,std::string>> inp_data = parse_input(argv[1]);
	std::map<std::string,std::string> params = std::get<1>(inp_data);
	System my_sys(std::get<0>(inp_data),std::get<1>(inp_data));
	my_sys.compute_cap_matrix();
	my_sys.CORRELATED_CAP_MAT.raw_print();
	auto t_end = std::chrono::high_resolution_clock::now();
	std::cout << "Wall time:" << std::chrono::duration<double, std::milli>(t_end-t_start).count() << std::endl;

	/*
// MOLCAS //////
	size_t nstates = 5;
	System sys1("example/N2.xyz","example/bas1.bas");
	arma::mat Smat(sys1.bs.num_carts(),sys1.bs.num_carts());
	compute_analytical_overlap(sys1.bs,Smat);
	uniform_cart_norm(Smat,sys1.bs);
	arma::mat spherical_ints(sys1.bs.Nbasis,sys1.bs.Nbasis);
	cart2spherical(Smat,spherical_ints,sys1.bs);
	to_molcas_ordering(spherical_ints, sys1.bs,sys1.atoms);
	auto dms = readRassiHDF5("example/caspt2_test.rassi.h5");
	arma::mat overlap_mat;
	overlap_mat.load(arma::hdf5_name("example/caspt2_test.rassi.h5", "AO_OVERLAP_MATRIX"));
	std::cout << std::endl;
	overlap_mat.reshape(sqrt(overlap_mat.n_cols),sqrt(overlap_mat.n_cols));
	std::cout << "Comparing matrix elements:" << std::endl;
	bool conflicts = false;
	for (size_t i=0;i<overlap_mat.n_rows;i++)
	{
		for(size_t j=0;j<overlap_mat.n_cols;j++)
		{
			if (abs(overlap_mat(i,j)-spherical_ints(i,j))>1E-5)
			{
				std::cout << "Conflict at:" << i << "," << j << std::endl;
				std::cout << "Molcas says:" << overlap_mat(i,j) << std::endl;
				std::cout << "Libcap says:" << spherical_ints(i,j) << std::endl;
				std::cout << abs(overlap_mat(i,j)-spherical_ints(i,j)) << std::endl;
				conflicts = true;
			}
		}
	}
	if (!conflicts)
		std::cout << "Overlap matrices match, good to go!" << std::endl;
	arma::mat cap_mat(sys1.bs.num_carts(),sys1.bs.num_carts());
	cap_mat.zeros();
	auto start = std::chrono::high_resolution_clock::now();
	compute_numerical_overlap_mat(cap_mat, sys1.bs, sys1.atoms);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	std::cout << "Time elapsed:" << duration.count() << std::endl;
	std::cout << "Finished calculating, now transforming..." << std::endl;
	uniform_cart_norm(cap_mat,sys1.bs);
	arma::mat cap_spherical(sys1.bs.Nbasis,sys1.bs.Nbasis);
	cart2spherical(cap_mat,cap_spherical,sys1.bs);
	cap_spherical.print();
	to_molcas_ordering(cap_spherical, sys1.bs,sys1.atoms);
	std::cout << "Finished the cap matrix." << std::endl;
	arma::mat EOMCAP(nstates,nstates);
	EOMCAP.zeros();
	for (size_t row_idx=0;row_idx<EOMCAP.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<EOMCAP.n_cols;col_idx++)
		{
			std::cout << "TDM for state:" << row_idx << "," << col_idx << std::endl;
			dms[row_idx][col_idx].print();
			EOMCAP(row_idx,col_idx) = -1*arma::trace(dms[row_idx][col_idx]*cap_spherical)-
									  arma::trace(dms[row_idx][col_idx]*cap_spherical);
		}
	}
	std::cout << std::endl;
	std::cout << std::fixed << std::setprecision(10);
	EOMCAP.raw_print(std::cout);
	*/
	return 0;
}

