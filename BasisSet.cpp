#include <iostream>
#include <fstream>
#include <vector>
#include "BasisFunction.h"
#include "BasisSet.h"
#include <math.h>
#include "Atom.h"
#include <sstream>
#include <algorithm>
#include <string>
# define M_PIl          3.141592653589793238462643383279502884
using namespace std;


BasisSet::BasisSet(std::string xyz_name,std::string basis_name)
{
	//carts = purecart;
	name = "user-specified";
	Nbasis = generateBasis(xyz_name,basis_name);
}

int BasisSet::generateBasis(std::string xyz_name,std::string basis_name)
{
	std::vector<Atom> geometry = read_xyz(xyz_name);
	return 0;
}

std::vector<Atom> BasisSet::read_xyz(std::string xyz_name)
{
	ifstream xyzfile;
	xyzfile.open(xyz_name);
	//read in number of atoms
	size_t natom;
	xyzfile >> natom;
	std::string rest_of_line;
	std::getline(xyzfile, rest_of_line);
	//read in comment
	std::string comment;
	std::getline(xyzfile, comment);
	//now the coordinates
	std::vector<Atom> atoms;
	for (size_t i=0;i<natom;i++)
	{
	    // read line
	    std::string linestr;
	    std::getline(xyzfile, linestr);
	    std::istringstream iss(linestr);
	    std::string element_symbol;
	    double x, y, z;
	    iss >> element_symbol >> x >> y >> z;
	    transform(element_symbol.begin(),element_symbol.end(),element_symbol.begin(),::toupper);
	    std::cout << element_symbol << std::endl;
	    atoms.push_back(Atom(element_symbol,x,y,z));
	}
	return atoms;


}
