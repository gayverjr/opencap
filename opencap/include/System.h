 /*! \file System.h
     \brief Class which holds molecular geometry and basis set data.
 */
#include <pybind11/pybind11.h>
#include "Atom.h"
#include "BasisSet.h"
#include <map>
#include <iostream>
#include <pybind11/numpy.h>
#include <eigen3/Eigen/Dense>

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
	System(){python=false;};
	/** Construct from python
	 */
	System(py::dict dict);
	/** Constructor from geometry and parameters
	 */
	System(std::vector<Atom> geometry,std::map<std::string, std::string> params);
	/**
	 * Constructor from output file
	 */
	System(std::string filename,std::string file_type);
	/** Overlap matrix in AO basis
	 */
	Eigen::MatrixXd OVERLAP_MAT;
	/** Sets geometry from python
	 */
	std::vector<Atom> parse_geometry_string(std::string geometry_string);
	/** Returns overlap matrix
	 */
	Eigen::MatrixXd get_overlap_mat();
	/** Compares computed overlap matrix with supplied overlap matrix
	 *  We also try to re-normalize, see renormalize_overlap
	 */
	bool check_overlap_mat(Eigen::MatrixXd smat, std::string ordering, std::string basis_file="");
	/** Set to true when constructed from the python interpreter, important for printing
	 */
	bool python;
	/** Returns ordered string rep of basis function IDs
	 */
	std::string get_basis_ids();


private:
	/** Converts atomic coordinates from angstrom to bohr units
	 */
	void convert_ang_to_bohr();
	/** Checks that parameters are sufficient to construct system,
	 * and sets unspecified fields to their defaults.
	 * \return Updated parameters map.
	 */
	void verify_system();
	void renormalize_overlap(Eigen::MatrixXd smat);
};



#endif /* SYSTEM_H_ */
