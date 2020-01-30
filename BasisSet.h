#include <vector>
#include "BasisFunction.h"
#include "Atom.h"
#include <map>
#include "Shell.h"
#pragma once
using namespace std;

class BasisSet {
public:
	std::string name;
	int Nbasis;
	//std::vector<int> carts;
	std::vector<BasisFunction> basis;

public:
	BasisSet(std::string xyz_name, std::string basis_name);
private:
	int generateBasisFunctions(std::string xyz_name,std::string basis_name);
	std::vector<Atom> read_xyz(std::string xyz_name);
	std::map<string,std::vector<Shell>> readBasis(std::string basis_name);

};
