 /*! \file BasisSet.h
     \brief Class containing data and helper functions pertaining to ab initio basis set.
 */
#include <vector>
#include "Atom.h"
#include <map>
#include "Shell.h"
#include <string>
#pragma once
using namespace std;


/*! \brief Contains data and helpful functions pertaining to the ab initio basis set.
 *
 *	The BasisSet class is essentially a wrapper around the vector of Shell objects which
 *	comprise the ab initio basis set. It contains a number of helpful functions such as computing the
 *	number of cartesian basis functions, computing the number of basis after rotation into
 *	spherical harmonic basis, maximum angular momentum etc.
 *
 *	Constructing a BasisSet object requires a geometry, and basis set parameters, including the basis set file
 *	containing primitive exponents and coefficients. Currently, only basis set files in Psi4 format are supported.
 *
 */

class BasisSet {
public:
	/** Number of gaussian shells specified in the ab initio basis set.
	 */
	size_t Nshells;
    /** List of Shell objects comprising the ab initio basis set.
      */
	std::vector<Shell> basis;
	   /** Number of basis functions after rotation in spherical harmonic basis.
	    */
	size_t Nbasis;

public:
    /** Default constructor.
      *
	  * For current version of OpenCAP, does nothing.
      */
	BasisSet();
    /** Default constructor.
      *
	  * Constructs BasisSet using geometry and dictionary of parameters which define ab initio basis set.
	  * \param geometry: molecular geometry
	  * \param parameters: parameters map
      */
	BasisSet(std::vector<Atom> geometry,std::map<std::string, std::string> parameters);
    /** Returns number of Cartesian basis functions in basis set.
      */
	size_t num_carts();
    /** Returns largest angular momentum of basis functions in the basis set.
      */
	int max_L();
    /** Returns smallest exponent of all primitives in the basis set centered at a particular Atom
     * for each orbital angular momentum.
     * \param atm: Atom of interest
     * \return Vector of smallest exponents belonging to a particular atom
      */
	std::vector<double> alpha_min(Atom atm);
    /** Returns largest exponent of all primitives in the basis set centered at a particular Atom.
     * \param atm: Atom of interest
     * \return Largest exponent in the basis set belonging to a particular atom
      */
	double alpha_max(Atom atm);
private:
    /** Calculates number of basis functions based on the Shells stored in the basis vector.
      */
	size_t calc_basis_size();
    /** Returns a list of Shells centered on a particular Atom.
     * \param atm: Atom of interest
     * \return Vector containing all Shell objects centered on the atom
      */
	std::vector<Shell> shells_on_center(Atom atm);
    /** Populates basis class member with Shells read in from the ab initio basis set. The coordinates of the Shell
     * are updated to the Atom it belongs to before being stored in the basis vector.
     * \param geometry: molecular geometry
     * \param all_shells: mapping of shells to atom types
      */
	void build_basis_set(std::vector<Atom> geometry,map<string,std::vector<Shell>> all_shells);
};
