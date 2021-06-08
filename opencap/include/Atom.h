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

/*! \file Atom.h
     \brief Classes and structs for storing atomic data.
 */

#include <array>
#include <string>
#include <vector>

#pragma once
using namespace std;

/*! \brief Object for storing element data
 *
 * Contains atomic name, symbol, and nuclear charge.
 *
 */
struct element
{
	/** Nuclear charge
	  */
    int Z;
    /** Element name
      */
    std::string name;
    /** Atomic symbol
      */
    std::string symbol;
    /**
     * \param new_Z: atomic number
     * \param new_name: element name
     * \param new_symbol: element symbol
     */
    element(int new_Z, std::string new_name, std::string new_symbol)
    {Z=new_Z; name = new_name; symbol=new_symbol;};
    /* Returns element name
     * \return name: Element name
     */
    std::string &getName() {return name;};
};


/*! \brief Object for storing atomic data
 *
 * Contains atomic name, symbol, nuclear charge, and coordinates.
 *
 */
class Atom {
public:
    /** Atomic symbol
      */
	std::string symbol;
    /** Element name
      */
	std::string name;
    /** Nuclear charge
      */
	int Z;
    /** Atomic coordinates, assumed to be in bohr units by the rest of the program
      */
	std::array<double,3> coords;

public:
    /** Converts coordinates from angstrom to bohr units.
      */
	void ang_to_bohr();
    /** Constructs Atom object from atomic symbol and coordinates.
     * \param symbol: element symbol
     * \param x: x coordinate
     * \param y: y coordinate
     * \param z: z coordinate
      */
	Atom(string symbol,double x,double y, double z);
    /** Constructs Atom object from atomic number and coordinates.
     * \param atom_num: atomic number
     * \param x: x coordinate
     * \param y: y coordinate
     * \param z: z coordinate
      */
	Atom(int atom_num,double x,double y, double z);
    /** Returns element struct corresponding to atomic symbol
     * \param symbol: atomic symbol
     * \return element: Element struct corresponding to atomic symbol
      */
	element get_element_info(std::string symbol);
    /** Returns element struct corresponding to atomic number
     * \param Z: atomic number
     * \return element: Element struct corresponding to atomic number
      */
	element get_element_info(int Z);
};
