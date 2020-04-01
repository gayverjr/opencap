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
#include "readMolcasHDF5.h"
#include "molcas_transform.h"
#include <chrono>
#include <iomanip>
#include <cmath>
#include <limits>
using namespace std;

int main(int argc, char **argv)
{
	// QCHEM ////
	size_t nstates = 6;
	System sys1("example/N2.xyz","example/big_bas.bas");
	std::cout << "Number of basis functions:" << sys1.bs.Nbasis << std::endl;
	//overlap matrix, for testing purposes
	arma::mat Smat(sys1.bs.num_carts(),sys1.bs.num_carts());
	compute_analytical_overlap(sys1.bs,Smat);
	uniform_cart_norm(Smat,sys1.bs);
	arma::mat spherical_ints(sys1.bs.Nbasis,sys1.bs.Nbasis);
	cart2spherical(Smat,spherical_ints,sys1.bs);
	to_molden_ordering(spherical_ints, sys1.bs);
	auto qchem_smat = qchem_read_overlap("example/re3.fchk",sys1.bs.Nbasis);
	std::cout << std::endl;
	std::cout << "Checking overlap matrix:" << std::endl;
	bool conflicts = false;
	for (size_t i=0;i<qchem_smat.n_rows;i++)
	{
		for(size_t j=0;j<qchem_smat.n_cols;j++)
		{
			if (abs(qchem_smat(i,j)-spherical_ints(i,j))>1E-5)
			{
				std::cout << "Conflict at:" << i << "," << j << std::endl;
				std::cout << "qchem says:" << qchem_smat(i,j) << std::endl;
				std::cout << "Libcap says:" << spherical_ints(i,j) << std::endl;
				std::cout << abs(qchem_smat(i,j)-spherical_ints(i,j)) << std::endl;
				conflicts = true;
			}
		}
	}
	if (conflicts)
		std::cout << "Should probably throw an error here." << std::endl;
	else
		std::cout << "Overlap matrices match, we're good to go!" << std::endl;
	std::cout << "Reading in the density matrices..." << std::endl;
	auto dms = qchem_read_in_dms("example/re3.fchk",nstates,sys1.bs.Nbasis);
	std::cout << "Finished reading in the density matrices..." << std::endl;
	//cap matrix stuff for later
	std::cout << "Starting cap matrix evaluation in AO basis." << std::endl;
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
	std::cout << "CAP matrix in AO basis" << std::endl;
	cap_spherical.print();
	to_molden_ordering(cap_spherical, sys1.bs);
	std::cout << "Finished the cap matrix." << std::endl;
	std::cout << std::endl;
	//read in DMs
	arma::mat EOMCAP(nstates,nstates);
	EOMCAP.zeros();
    std::cout << std::fixed << std::setprecision(10);
	for (size_t row_idx=0;row_idx<EOMCAP.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<EOMCAP.n_cols;col_idx++)
		{
			EOMCAP(row_idx,col_idx) = arma::trace(dms[0][row_idx][col_idx]*cap_spherical)+
									  arma::trace(dms[1][row_idx][col_idx]*cap_spherical);
			EOMCAP(row_idx,col_idx) = -1.0* EOMCAP(row_idx,col_idx);

		}
	}
	std::cout << "Final cap matrix" << std::endl;
	EOMCAP.raw_print(std::cout);




// MOLCAS //////
	/*
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

