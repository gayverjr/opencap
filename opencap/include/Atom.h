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
    /**
     * \param new_Z: atomic number
     * \param new_name: element name
     * \param new_symbol: element symbol
     */
    element(int new_Z, std::string new_name, std::string new_symbol)
    {Z=new_Z; name = new_name; symbol=new_symbol;};
    /* Returns element name
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
    /** Returns element struct corresponding to atomic symbol.
      */
	element get_element_info(std::string symbol);
    /** Returns element struct corresponding to atomic number.
      */
	element get_element_info(int Z);
};
