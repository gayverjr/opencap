#include <vector>
#include "Atom.h"
#include <map>
#include "Shell.h"
#pragma once
using namespace std;

class BasisSet {
public:
	std::string name;
	size_t Nshells;
	std::vector<Shell> basis;
	size_t Nbasis;

public:
	BasisSet(std::string xyz_name, std::string basis_name);
	size_t num_carts();
private:
	size_t generateBasisFunctions(std::string xyz_name,std::string basis_name);
	size_t calc_basis_size();
	std::vector<Atom> read_xyz(std::string xyz_name);
	std::map<string,std::vector<Shell>> readBasis(std::string basis_name);

};
