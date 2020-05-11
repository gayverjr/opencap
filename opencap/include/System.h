 /*! \file System.h
     \brief Class which holds all of the data pertaining to the physical system.
 */
#include "Atom.h"
#include "BasisSet.h"
#include "CAP.h"
#include <map>
#include <armadillo>
#ifndef SYSTEM_H_
#define SYSTEM_H_


/*! \brief Class which holds all data pertaining to the physical system.
 *
 * Currently, this class is capable of computing the %CAP in AO basis, and with data
 * read in from the electronic structure package, in correlated many electron basis.
 *
 */

class System
{
public:
	/** Atoms in the system
	 */
	std::vector<Atom> atoms;
	/** Map containing the parameters defined in the input
	 */
	std::map<std::string, std::string> parameters;
	/** Ab initio basis set
	 */
	BasisSet bs;
	/** Default constructor, does nothing
	 */
	System();
	/** Constructor from geometry and paramters
	 */
	System(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	/** %CAP matrix in AO basis
	 */
	arma::mat AO_CAP_MAT;
	/** %CAP matrix in correlated many electron basis
	 */
	arma::mat CORRELATED_CAP_MAT;
	/** Zeroth order Hamiltonian. Dimension is (nstates,nstates)
	 */
	arma::mat ZERO_ORDER_H;
	/** Overlap matrix in AO basis
	 */
	arma::mat OVERLAP_MAT;
	/** Transition density matrices in AO basis, alpha densities
	 */
	std::vector<std::vector<arma::mat>> alpha_dms;
	/** Transition density matrices in AO basis, beta densities
	 */
	std::vector<std::vector<arma::mat>> beta_dms;
	/** Number of states
	 */
	size_t nstates;
	/** Computes %CAP in AO basis, and then correlated basis, saves to respective class members
	 */
	void compute_cap_matrix();

private:
	/** Compares computed overlap matrix to that read in from the electronic structure package
	 */
	bool check_overlap_matrix();
	/** Converts atomic coordinates from angstrom to bohr units
	 */
	void convert_ang_to_bohr();
	/** Checks %CAP section of input to ensure that all required keywords have valid values
	 */
	bool verify_cap_parameters(std::string key);
	/** Checks %System section of input to ensure that all required keywords have valid values
	 */
	bool check_parameters();
	/** Reads in TDMs from electronic structure package
	 */
	bool read_in_dms();
	/** Computes %CAP in correlated many electron basis
	 */
	void compute_cap_correlated_basis();
	/** Re-orders %CAP matrix in AO basis to match electronic structure package ordering
	 */
	void reorder_cap();
	/** Verifies choice of electronic structure method
	 */
	bool verify_method(std::string key);
	/** Reads in zeroth order Hamiltonian from electronic structure package
	 */
	bool read_in_zero_order_H();
	/** Reads in zeroth order Hamiltonian from file
	 */
	arma::mat read_h0_file();
};



#endif /* SYSTEM_H_ */
