/*Copyright (c) 2021 James Gayvert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*! \file System.h
     \brief Class which holds molecular geometry and basis set data.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <Eigen/Dense>
#include <iostream>
#include <map>

#include "BasisSet.h"


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
