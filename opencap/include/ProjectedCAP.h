 /*! \file ProjectedCAP.h
     \brief Class which handles projected %CAP calculations.
 */
#include <pybind11/pybind11.h>
#include "opencap_exception.h"
#include "System.h"
#include "CAP.h"
#include <armadillo>
#include <pybind11/numpy.h>

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
	arma::mat AO_CAP_MAT;
	/** %CAP matrix in correlated many electron basis
	 */
	arma::mat CORRELATED_CAP_MAT;
	bool python;
	/** Zeroth order Hamiltonian. Dimension is (nstates,nstates)
	 */
	arma::mat ZERO_ORDER_H;
	/** Transition density matrices in AO basis, alpha densities
	 */
	std::vector<std::vector<arma::mat>> alpha_dms;
	/** Transition density matrices in AO basis, beta densities
	 */
	std::vector<std::vector<arma::mat>> beta_dms;
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
	void compute_cap_matrix();
	/** Checks that electronic structure method and package is supported, and that necessary keywords are present.
	 */
	void verify_method();
	/** Executes Projected CAP method.
	 */
	void run();
	Projected_CAP(System my_sys,size_t num_states,std::string gto_ordering);
	py::array get_AO_CAP();
	py::array get_CAP_mat();
	py::array get_H();
	void add_tdm(py::array_t<double> & alpha_density,
			py::array_t<double> & beta_density,size_t row_idx, size_t col_idx);
	void add_tdm(py::array_t<double> & tdm, size_t row_idx, size_t col_idx);
	void set_h0(py::array_t<double> &h0);
	void read_electronic_structure_data(py::dict dict);
	void set_cap_params(py::dict dict);

private:
	/** Reads in TDMs from electronic structure package
	 */
	void read_in_dms();
	/** Computes %CAP in correlated many electron basis
	 */
	void compute_cap_correlated_basis();
	/** Re-orders %CAP matrix in AO basis to match electronic structure package ordering
	 */
	void reorder_cap();
	/** Reads in zeroth order Hamiltonian from electronic structure package
	 */
	void read_in_zero_order_H();
	/** Reads in zeroth order Hamiltonian from file
	 */
	arma::mat read_h0_file();
	/** Compares computed overlap matrix to that read in from the electronic structure package
	 */
	void check_overlap_matrix();
};

#endif /* INCLUDE_PROJECTEDCAP_H_ */
