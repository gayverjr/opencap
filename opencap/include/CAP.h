 /*! \file CAP.h
     \brief Class which handles computing the %CAP matrix.
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "opencap_exception.h"
#include "System.h"
#include "AOCAP.h"
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#ifndef INCLUDE_CAP_H_
#define INCLUDE_CAP_H_

/*! \brief Class which handles computing the %CAP matrix.
 *
 */
class CAP
{
public:
	/** Overlap matrix parsed from file, can be empty when constructed from python
	 */
	Eigen::MatrixXd OVERLAP_MAT;
	/** Geometry and basis set
	 */
	System system;
	/** %CAP matrix in AO basis
	 */
	Eigen::MatrixXd AO_CAP_MAT;
	/** %CAP matrix in correlated many electron basis
	 */
	Eigen::MatrixXd CAP_MAT;
	/** Set to true when constructed from the python interpreter, important for printing
	 */
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
	/** Constructs %CAP object from %System object.
	 *  \param my_sys: System object
	 *  \param params: Map of parameters
	 */
	CAP(System &my_sys,std::map<std::string, std::string> params);
	/** Computes %CAP in AO basis
	 */
	void compute_ao_cap();
	/** Computes %CAP in state basis
	 */
	void compute_perturb_cap();
	/** Checks that electronic structure method and package is supported, and that necessary keywords are present.
	 */
	void verify_method(std::map<std::string,std::string> params);
	/** Executes Perturbative CAP method.
	 */
	void run();
	/** Constructor for Python.
	 *\param my_sys: System object
	 *\param num_states: number of states
	 *\param gto_ordering: Name of electronic structure package
	 *\param dict: Python dictionary containing parameters
	 */
	CAP(System my_sys,py::dict dict,size_t num_states,std::string gto_ordering);
	/** Returns CAP matrix in AO basis.
	 */
	Eigen::MatrixXd get_ao_cap();
	/** Returns CAP matrix in wave function basis.
	 */
	Eigen::MatrixXd get_perturb_cap();
	/** Returns zeroth order Hamiltonian.
	 */
	Eigen::MatrixXd get_H();
	/** Adds transition density matrices (alpha and beta) from state row_idx --> col_idx.
	 * \param alpha_density: TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param beta_density: TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param row_idx: initial state index
	 * \param col_idx: final state index
	 * \param ordering: order of GTOs
	 * \param basis_file: File containing basis set specification. Required for OpenMolcas.
	 */
	void add_tdms(Eigen::MatrixXd & alpha_density,
			Eigen::MatrixXd & beta_density,size_t row_idx, size_t col_idx,
			std::string ordering, std::string basis_file="");
	/** Adds spin traced transition density matrix from state row_idx --> col_idx.
	 * \param tdm: Spin traced TDM in AO basis of dimension (NBasis,Nbasis)
	 * \param row_idx: initial state index
	 * \param col_idx: final state index
	 * \param ordering: order of GTOs
	 * \param basis_file: File containing basis set specification. Required for OpenMolcas.
	 */
	void add_tdm(Eigen::MatrixXd tdm,size_t row_idx, size_t col_idx,
			std::string ordering,std::string basis_file="");
	/** Reads in electronic structure data from file, from python.
	 * Valid keywords: method,qc_output,h0_file,rassi_h5,
			fchk_file,molcas_output.
	 * \param dict: Python dictionary containing keywords for reading in data from disk.
	 */
	void read_electronic_structure_data(py::dict dict);
	/** Verifies that required electronic structure data is present to perform calculation.
	 */
	void verify_data();
	/** Renormalizes the AO %CAP matrix using the supplied overlap matrix.
	 */
	void renormalize_cap(Eigen::MatrixXd smat, std::string ordering, std::string basis_file="");
	/** Renormalizes the AO %CAP matrix using the previously parsed overlap matrix.
	 */
	void renormalize();

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

#endif /* INCLUDE_CAP_H_ */
