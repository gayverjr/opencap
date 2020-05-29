 /*! \file System.h
     \brief Class which holds molecular geometry and basis set data.
 */
#include <pybind11/pybind11.h>
#include "Atom.h"
#include "BasisSet.h"
#include "CAP.h"
#include <map>
#include <armadillo>
#include <iostream>

#ifndef SYSTEM_H_
#define SYSTEM_H_
namespace py = pybind11;

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
	/** Construct from python
	 */
	System(py::dict dict);
	/** Constructor from geometry and parameters
	 */
	System(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	/** Overlap matrix in AO basis
	 */
	arma::mat OVERLAP_MAT;
	/** Sets geometry from python
	 */
	void set_geometry(std::string geometry_string);


private:
	/** Converts atomic coordinates from angstrom to bohr units
	 */
	void convert_ang_to_bohr();
	/** Checks that parameters are sufficient to construct system,
	 * and sets unspecified fields to their defaults.
	 * \return Updated parameters map.
	 */
	void verify_system_parameters(std::map<std::string, std::string> &params);
	/** Checks bohr_coordinates keyword.
	 */
	bool bohr_coords();
};



#endif /* SYSTEM_H_ */
