 /*! \file Atom.h
     \brief Classes and structs for storing atomic data.
 */
#include <vector>
#include <array>
#include <string>
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
    /** Returns element structure corresponding to atomic symbol.
      */
	element get_element_info(std::string symbol);
};
