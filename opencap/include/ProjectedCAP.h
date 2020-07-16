 /*! \file ProjectedCAP.h
     \brief Class which handles projected %CAP calculations.
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "opencap_exception.h"
#include "System.h"
#include "CAP.h"
#include <pybind11/eigen.h>

#ifndef INCLUDE_PROJECTEDCAP_H_
#define INCLUDE_PROJECTEDCAP_H_

/*! \brief Class which handles projected %CAP calculations.
 *
 */
class Projected_CAP
{
public:
	System system;
	/** %CAP matrix in AO basis
	 */
	Eigen::MatrixXd AO_CAP_MAT;
	/** %CAP matrix in correlated many electron basis
	 */
	Eigen::MatrixXd CORRELATED_CAP_MAT;
	bool python;
	/** Zeroth order Hamiltonian. Dimension is (nstates,nstates)
	 */
	Eigen::MatrixXd ZERO_ORDER_H;
	/** Transition density matrices in AO basis, alpha densities
	 */
	std::vector<std::vector<Eigen::MatrixXd>> alpha_dms;
	/** Transition density matrices in AO basis, beta densities
	 */
	std::vector<std::vector<Eigen::MatrixXd>> beta_dms;
	/** Number of states
	 */
	size_t nstates;
	/** Map containing the parameters defined in the input
	 */
	std::map<std::string, std::string> parameters;
	/** Constructs Projected %CAP object from %System object.
	 *  \param my_sys: System object
	 */
	Projected_CAP(System my_sys,std::map<std::string, std::string> params);
	/** Computes %CAP in AO basis, and then correlated basis, saves to respective class members
	 */
	void compute_ao_cap();
	/** Computes %CAP in correlated many electron basis
	 */
	void compute_projected_cap();
	/** Checks that electronic structure method and package is supported, and that necessary keywords are present.
	 */
	void verify_method(std::map<std::string,std::string> params);
	/** Executes Projected CAP method.
	 */
	void run();
	/** Constructor for Python.
	 *\param my_sys: System object
	 *\param num_states: number of states
	 *\param gto_ordering: Name of electronic structure package
	 */
	Projected_CAP(System my_sys,py::dict dict,size_t num_states,std::string gto_ordering);
	/** Returns CAP matrix in AO basis.
	 */
	Eigen::MatrixXd get_ao_cap();
	/** Returns CAP matrix in wave function basis.
	 */
	Eigen::MatrixXd get_projected_cap();
	/** Returns zeroth order Hamiltonian in wave function basis.
	 */
	Eigen::MatrixXd get_H();
	/** Adds transition density matrices (alpha and beta) from state row_idx --> col_idx.
	 * \param alpha_density: TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param beta_density: TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param row_idx: initial state index
	 * \param col_idx: final state index
	 */
	void add_tdms(Eigen::MatrixXd & alpha_density,
			Eigen::MatrixXd & beta_density,size_t row_idx, size_t col_idx,
			std::string ordering, std::string basis_file="");
	/** Adds spin traced transition density matrix from state row_idx --> col_idx.
	 * \param alpha_density: Spin traced TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param row_idx: initial state index
	 * \param col_idx: final state index
	 */
	void add_tdm(Eigen::MatrixXd tdm,size_t row_idx, size_t col_idx,
			std::string ordering,std::string basis_file="");
	/** Reads in electronic structure data from file, from python.
	 * Valid keywords: method,qc_output,h0_file,rassi_h5,
			fchk_file,molcas_output.
	 * \param dict: Python dictionary containing keywords for reading in data from disk.
	 * \param override: Set to true when overriding existing keywords
	 */
	void read_electronic_structure_data(py::dict dict);
	/** Verifies that required electronic structure data is present to perform calculation.
	 */
	void verify_data();

private:
	/** Reads in TDMs from electronic structure package
	 */
	void read_in_dms();
	/** Reads in zeroth order Hamiltonian from electronic structure package
	 */
	void read_in_zero_order_H();
	/** Reads in zeroth order Hamiltonian from file
	 */
	Eigen::MatrixXd read_h0_file();
	/** Compares computed overlap matrix to that read in from the electronic structure package
	 */
	void check_overlap_matrix();
};

#endif /* INCLUDE_PROJECTEDCAP_H_ */
