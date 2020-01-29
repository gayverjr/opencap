#include <vector>
#include "BasisFunction.h"
#include "Atom.h"
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
	int generateBasis(std::string xyz_name,std::string basis_name);
	std::vector<Atom> read_xyz(std::string xyz_name);

};
