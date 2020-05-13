 /*! \file System.h
     \brief Class which holds molecular geometry and basis set data.
 */
#include "Atom.h"
#include "BasisSet.h"
#include "CAP.h"
#include <map>
#include <armadillo>
#ifndef SYSTEM_H_
#define SYSTEM_H_


/*! \brief Class which holds molecular geometry and basis set data.
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
	/** Constructor from geometry and parameters
	 */
	System(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	/** Overlap matrix in AO basis
	 */
	arma::mat OVERLAP_MAT;

private:
	/** Converts atomic coordinates from angstrom to bohr units
	 */
	void convert_ang_to_bohr();
	/** Checks that parameters are sufficient to construct system.
	 */
	void verify_system_parameters(std::map<std::string, std::string> &params);
	/** Checks bohr_coordinates keyword.
	 */
	bool bohr_coords();
};



#endif /* SYSTEM_H_ */
