 /*! \file ProjectedCAP.h
     \brief Class which handles projected %CAP calculations.
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "opencap_exception.h"
#include "System.h"
#include "CAP.h"
#include <armadillo>

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
	void verify_method(std::map<std::string,std::string> params);
	/** Executes Projected CAP method.
	 */
	void run();
	/** Constructor for Python.
	 *\param my_sys: System object
	 *\param num_states: number of states
	 *\param gto_ordering: Name of electronic structure package
	 */
	Projected_CAP(System my_sys,size_t num_states,std::string gto_ordering);
	/** Returns CAP matrix in AO basis.
	 */
	py::array get_AO_CAP();
	/** Returns CAP matrix in wave function basis.
	 */
	py::array get_CAP_mat();
	/** Returns zeroth order Hamiltonian in wave function basis.
	 */
	py::array get_H();
	/** Adds transition density matrices (alpha and beta) from state row_idx --> col_idx.
	 * \param alpha_density: TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param beta_density: TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param row_idx: initial state index
	 * \param col_idx: final state index
	 */
	void add_tdm(py::array_t<double> & alpha_density,
			py::array_t<double> & beta_density,size_t row_idx, size_t col_idx);
	/** Adds spin traced transition density matrix from state row_idx --> col_idx.
	 * \param alpha_density: Spin traced TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param row_idx: initial state index
	 * \param col_idx: final state index
	 */
	void add_tdm(py::array_t<double> & tdm, size_t row_idx, size_t col_idx);
	/** Sets zeroth order Hamiltonian.
	 * \param h0: Zeroth order Hamiltonian of dimension (nstates,nstates)
	 */
	void set_h0(py::array_t<double> &h0);
	/** Reads in electronic structure data from file, from python.
	 * Valid keywords: method,qc_output,h0_file,rassi_h5,
			fchk_file,molcas_output.
	 * \param dict: Python dictionary containing keywords for reading in data from disk.
	 * \param override: Set to true when overriding existing keywords
	 */
	void read_electronic_structure_data(py::dict dict);
	/** Sets CAP parameters from python.
	 * Valid keywords: cap_type,cap_x,cap_y,cap_z,
	 * r_cut,radial_precision,angular_points.
	 * \param dict: Python dictionary containing keywords for setting CAP parameters.
	 * \param override: Set to true when overriding existing keywords
	 */
	void set_cap_params(py::dict dict);
	void verify_data();

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
