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
	System(){};
	/** Constructor from geometry and paramters
	 */
	System(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	/** Overlap matrix in AO basis
	 */
	arma::mat OVERLAP_MAT;

private:
	/** Converts atomic coordinates from angstrom to bohr units
	 */
	void convert_ang_to_bohr();
	void verify_system_parameters(std::map<std::string, std::string> &params);
};



#endif /* SYSTEM_H_ */
