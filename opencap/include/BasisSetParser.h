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

/*! \file BasisSetParser.h
     \brief Class for parsing basis set files.
 */

#include <map>
#include <string>
#include <vector>

#include "Shell.h"

#ifndef BASISSETPARSER_H_
#define BASISSETPARSER_H_

/*! \brief Parses the ab initio basis set from a basis set file.
 *
 *	The BasisSetParser class is used by the BasisSet class to read in an ab initio basis set. Basis sets augmented
 *	with extra sets of diffuse functions are typically used for CAP augmented calculations. We suggest downloading
 *	the basis set file from https://www.basissetexchange.org/ in Psi4 format, and then modifying the file to add
 *	the additional basis functions. Currently, the class only parses files in the Psi4 format.
 *
 *
 */

class BasisSetParser
{
public:
	/** String corresponding to which basis functions are cartesian/harmonic. For example,
	 * the string "spf" means that S,P, and F functions are Cartesian, while all other types are spherical
	 * harmonic.
	 */
	std::string cart_bf;
	/** String location of basis set file to be read
	 */
	std::string basis_file;
	BasisSetParser(std::map<std::string, std::string> parameters);
	/*! Read in basis set from file.
	    \return Mapping of atomic symbol to a list of shells corresponding to that atom in the basis set.
	*/
	std::map<std::string,std::vector<Shell>> read_basis();


private:
	/*! Read in basis set from file.
	    \return Mapping of atomic symbol to a list of shells corresponding to that atom in the basis set.
	*/
	std::map<std::string,std::vector<Shell>> read_basis_file();

};


#endif /* BASISSETPARSER_H_ */
