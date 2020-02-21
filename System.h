/*
 * System.h
 *
 *  Created on: Feb 18, 2020
 *      Author: JG
 */
#include "Atom.h"
#include "BasisSet.h"
#include "CAP.h"
#ifndef SYSTEM_H_
#define SYSTEM_H_

class System
{
public:
	int charge;
	int spin;
	std::vector<Atom> atoms;
	//std::string cart;
	BasisSet bs;
	//int nelectron;
	System(std::string xyz_name,std::string basis_name);

private:
	std::vector<Atom> read_xyz(std::string xyz_name);
	std::map<string,std::vector<Shell>> readBasis(std::string basis_name);

};



#endif /* SYSTEM_H_ */
