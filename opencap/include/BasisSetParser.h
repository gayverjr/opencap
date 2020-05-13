 /*! \file BasisSetParser.h
     \brief Class for parsing basis set files.
 */
#include <string>
#include <map>
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
	map<string,std::vector<Shell>> read_basis();


private:
	/*! Read in basis set from file.
	    \return Mapping of atomic symbol to a list of shells corresponding to that atom in the basis set.
	*/
	map<string,std::vector<Shell>> read_basis_file();

};


#endif /* BASISSETPARSER_H_ */
