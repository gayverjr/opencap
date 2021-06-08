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

/*! \file BasisSet.h
     \brief Classes and structs for storing basis set data.
 */
#pragma once

#include <map>
#include <string>
#include <vector>

#include "Shell.h"

using namespace std;


/*! \brief Struct for ID of %Shell
 *
 * Contains center (atom) index, principal quantum number n, angular momentum quantum number l. Angular momentum
 * quantum number l is negative when the shell is cartesian.
 *
 */
struct shell_id
{
	/** Index of atom this shell belongs to
	 */
	size_t ctr;
	/* Shell number
	 */
	size_t shell_num;
	/** Angular momentum quantum number. Positive if pure, negative if Cartesian
	 */
	int l;
	shell_id(size_t center,size_t shell_n, int l_num)
	{ctr=center;shell_num=shell_n;l=l_num;};
	bool operator==(const shell_id& other);
};

/*! \brief Struct for ID of basis function
 *
 * Contains center (atom) index, principal quantum number n, angular momentum quantum number l, magnetic quantum number m.
 * Angular momentum quantum number l is negative when the shell is cartesian.
 */
struct bf_id
{
	/** Index of atom this shell belongs to
	 */
	size_t ctr;
	/** Shell number
	 */
	size_t shell_num;
	/** Angular momentum quantum number. Positive if pure, negative if Cartesian
	 */
	int l;
	/** Magnetic quantum number. For cartesian basis functions, it goes 0 to (l + 1) * (l + 2) / 2
	 */
	int m;
	bf_id(shell_id id,int new_m)
	{ctr=id.ctr;shell_num=id.shell_num;l=id.l;m=new_m;};
	bool operator==(const bf_id& other);
};


/*! \brief Object for storing basis set information.
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
	/* Ordered list of shell_ids.
	 */
	std::vector<shell_id> shell_ids;
	/* Ordered list of bf_ids, specifies the ordering of the basis set.
	 */
	std::vector<bf_id> bf_ids;
	/* List of centers in basis set.
	 */
	std::vector<std::array<double,3>> centers;

public:
    /** Default constructor.
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
     * \return num_carts: Number of basis functions in cartesian representation
      */
	size_t num_carts();
    /** Returns largest angular momentum of basis functions in the basis set.
     * \return max_L: largest angular momentum in the basis set
      */
	int max_L();
    /** Returns smallest exponent of all primitives in the basis set centered at a particular Atom
     * for each orbital angular momentum.
     * \param atm: Atom of interest
     * \return Vector of smallest exponents belonging to a particular atom
      */
	std::vector<double> alpha_min(Atom &atm);
    /** Returns largest exponent of all primitives in the basis set centered at a particular Atom.
     * \param atm: Atom of interest
     * \return Largest exponent in the basis set belonging to a particular atom
      */
	double alpha_max(Atom &atm);
	/** Adds a new shell to the basis set and updates the data structures.
	 * \param new_shell: New Shell object to be added to the basis set.
	 */
	void add_shell(Shell &new_shell);
	/** Gets index in shell_ids list of given shell_id.
	 * \param id: Shell id
	 * \return  Index of shell in the basis set
	 */
	long get_index_of_shell_id(shell_id id);
	/** Normalizes all of the shells in the basis set.
	 */
	void normalize();
private:
    /** Calculates number of basis functions based on the Shells stored in the basis vector.
      */
	size_t calc_basis_size();
    /** Returns a list of Shells centered on a particular Atom.
     * \param atm: Atom of interest
     * \return Vector containing all Shell objects centered on the atom
      */
	std::vector<Shell> shells_on_center(Atom &atm);
    /** Populates basis class member with Shells read in from the ab initio basis set. The coordinates of the Shell
     * are updated to the Atom it belongs to before being stored in the basis vector.
     * \param geometry: molecular geometry
     * \param all_shells: mapping of shells to atom types
      */
	void build_basis_set(std::vector<Atom> geometry,map<string,std::vector<Shell>> all_shells);
};
