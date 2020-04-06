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
	BasisSet();
	BasisSet(std::vector<Atom> geometry,std::map<std::string, std::string> parameters);
	size_t num_carts();
	size_t max_L();
	std::vector<double> alpha_min(Atom atm);
	double alpha_max(Atom atm);
private:
	size_t generateBasisFunctions(std::vector<Atom> geometry,map<string,std::vector<Shell>> all_shells);
	size_t calc_basis_size();
	std::vector<Shell> shells_on_center(Atom atm);
	void build_basis_set(std::vector<Atom> geometry,map<string,std::vector<Shell>> all_shells);
};
