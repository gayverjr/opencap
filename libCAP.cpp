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
using namespace std;

int main(int argc, char **argv)
{
	//small system example
	std::cout << "Small system" << std::endl;
	size_t nstates = 3;
	System sys1("example/N2.xyz","example/bas1.bas");
	std::cout << "Number of basis functions:" << sys1.bs.Nbasis << std::endl;
	//overlap matrix, for testing purposes
	arma::mat Smat(sys1.bs.num_carts(),sys1.bs.num_carts());
	compute_analytical_overlap(sys1.bs,Smat);
	uniform_cart_norm(Smat,sys1.bs);
	arma::mat spherical_ints(sys1.bs.Nbasis,sys1.bs.Nbasis);
	cart2spherical(Smat,spherical_ints,sys1.bs);
	to_molden_ordering(spherical_ints, sys1.bs);
	auto qchem_smat = qchem_read_overlap("example/cc_test.fchk",sys1.bs.Nbasis);
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
	to_molden_ordering(cap_spherical, sys1.bs);
	std::cout << "Finished the cap matrix." << std::endl;
	std::cout << std::endl;
	//read in DMs
	std::cout << "Reading in the density matrices..." << std::endl;
	auto dms = qchem_read_in_dms("example/cc_test.fchk",nstates,sys1.bs.Nbasis);
	std::cout << "Finished reading in the density matrices..." << std::endl;
	std::cout << "Printing the CAP matrix!" << std::endl;
	arma::mat EOMCAP(nstates,nstates);
	for (size_t row_idx=0;row_idx<EOMCAP.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<EOMCAP.n_cols;col_idx++)
		{
			EOMCAP(row_idx,col_idx) = arma::trace(dms[0][row_idx][col_idx]*cap_spherical)+
									  arma::trace(dms[1][row_idx][col_idx]*cap_spherical);
		}
	}
	EOMCAP.print();



	//Big system example
	std::cout << "Big system" << std::endl;
	System my_system("example/N2.xyz","example/bas2.bas");
	std::cout << "Number of basis functions:" << my_system.bs.Nbasis << std::endl;
	//overlap matrix, for testing purposes
	arma::mat Smat2(my_system.bs.num_carts(),my_system.bs.num_carts());
	compute_analytical_overlap(my_system.bs,Smat2);
	uniform_cart_norm(Smat2,my_system.bs);
	arma::mat spherical_ints2(my_system.bs.Nbasis,my_system.bs.Nbasis);
	cart2spherical(Smat2,spherical_ints2,my_system.bs);
	to_molden_ordering(spherical_ints2, my_system.bs);
	auto qchem_smat2 = qchem_read_overlap("example/real.fchk",my_system.bs.Nbasis);
	std::cout << std::endl;
	std::cout << "Checking overlap matrix:" << std::endl;
	conflicts = false;
	for (size_t i=0;i<qchem_smat2.n_rows;i++)
	{
		for(size_t j=0;j<qchem_smat2.n_cols;j++)
		{
			if (abs(qchem_smat2(i,j)-spherical_ints2(i,j))>1E-5)
			{
				std::cout << "Conflict at:" << i << "," << j << std::endl;
				std::cout << "qchem says:" << qchem_smat2(i,j) << std::endl;
				std::cout << "Libcap says:" << spherical_ints2(i,j) << std::endl;
				std::cout << abs(qchem_smat2(i,j)-spherical_ints2(i,j)) << std::endl;
				conflicts = true;
			}
		}
	}
	if (conflicts)
		std::cout << "Should probably throw an error here." << std::endl;
	else
		std::cout << "Overlap matrices match, we're good to go!" << std::endl;

	//cap matrix stuff for later
	std::cout << "Starting cap matrix evaluation in AO basis." << std::endl;
	arma::mat cap_mat2(my_system.bs.num_carts(),my_system.bs.num_carts());
	cap_mat2.zeros();
	start = std::chrono::high_resolution_clock::now();
	compute_numerical_overlap_mat(cap_mat2, my_system.bs, my_system.atoms);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	std::cout << "Time elapsed:" << duration.count() << std::endl;
	std::cout << "Finished calculating, now transforming..." << std::endl;
	//end of q-chem stuff



	/* MOLCAS stuff, still in testing...
	//now hdf5 and molcas test
	auto dms2 = readRassiHDF5("ras66.rassi.h5");
    //to_molcas_ordering(spherical_ints, my_system.bs,my_system.atoms);
	//molcas overlap matrix for testing purposes
	arma::cube rassi_data;
	rassi_data.load(arma::hdf5_name("ras66.rassi.h5", "SFS_TRANSITION_DENSITIES"));
	arma::mat overlap_mat;
	overlap_mat.load(arma::hdf5_name("ras66.rassi.h5", "AO_OVERLAP_MATRIX"));
	std::cout << std::endl;
	overlap_mat.reshape(sqrt(overlap_mat.n_cols),sqrt(overlap_mat.n_cols));
	std::cout << "Comparing matrix elements:" << std::endl;
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
			}
		}
	}
	*/
	return 0;
}

