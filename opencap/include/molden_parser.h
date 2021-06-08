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

/*! \file molden_parser.h
     \brief Functions pertaining to reading the geometry and basis set from Molden.
     The format is detailed here: http://cheminf.cmbi.ru.nl/molden/molden_format.html
 */


#ifndef INCLUDE_MOLDEN_PARSER_H_
#define INCLUDE_MOLDEN_PARSER_H_

#include <string>

#include "BasisSet.h"


/*! Reads in geometry from molden file.
 * \param filename: file location of molden file
 * \return Geometry specified by molden file
*/
std::vector<Atom> read_geometry_from_molden(std::string filename);
/*! Reads in basis set from molden file.
 *  \param filename: file location of molden file
 *  \param atoms: Ordered list of atoms in the system
 *  \return BasisSet specified by molden file
*/
BasisSet read_basis_from_molden(std::string filename,std::vector<Atom> atoms);




#endif /* INCLUDE_MOLDEN_PARSER_H_ */
