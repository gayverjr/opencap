 /*! \file molden_parser.h
     \brief Functions pertaining to reading the geometry and basis set from Molden.
     The format is detailed here: http://cheminf.cmbi.ru.nl/molden/molden_format.html
 */
#include "Atom.h"
#include "BasisSet.h"


#ifndef INCLUDE_MOLDEN_PARSER_H_
#define INCLUDE_MOLDEN_PARSER_H_


std::vector<Atom> read_geometry_from_molden(std::string filename);
BasisSet read_basis_from_molden(std::string filename,std::vector<Atom> atoms);




#endif /* INCLUDE_MOLDEN_PARSER_H_ */
