/*
 * BasisSetParser.h
 *
 *  Created on: Apr 6, 2020
 *      Author: JG
 */
#include <string>
#include <map>
#include "Shell.h"

#ifndef BASISSETPARSER_H_
#define BASISSETPARSER_H_

class BasisSetParser
{
public:
	std::string cart_bf;
	std::string basis_file;
	std::string basis_name;
	BasisSetParser(std::map<std::string, std::string> parameters);
	map<string,std::vector<Shell>> read_basis();


private:
	map<string,std::vector<Shell>> read_basis_file(std::string basis_file);

};


#endif /* BASISSETPARSER_H_ */
