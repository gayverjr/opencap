#include <iostream>
#include <fstream>
#include <vector>
#include "BasisSet.h"
#include "Atom.h"
#include <math.h>
#include <sstream>
#include <algorithm>
#include <string>
#include <map>
#include <list>
#include "Shell.h"
#include <locale>
#include <iterator>
#include <cassert>
# define M_PIl          3.141592653589793238462643383279502884
using namespace std;

BasisSet::BasisSet(std::string xyz_name,std::string basis_name)
{
	//carts = purecart;
	name = "user-specified";
	Nshells = generateBasisFunctions(xyz_name,basis_name);
	Nbasis = calc_basis_size();
	std::cout << "Shells:" << Nshells << std::endl;
	std::cout << "Basis functions:" << Nbasis << std::endl;
	std::cout << "Number of cartesians:" << num_carts() << std::endl;

}

size_t BasisSet::calc_basis_size()
{
	size_t num_functions = 0;
	for(size_t i=0;i<Nshells;i++)
	{
		num_functions+=basis[i].num_bf;
	}
	return num_functions;
}

size_t BasisSet::num_carts()
{
	size_t num_carts = 0;
	for(auto&shell:basis)
	{
		num_carts += shell.num_carts();
	}
	return num_carts;
}

size_t BasisSet::generateBasisFunctions(std::string xyz_name,std::string basis_name)
{
	std::vector<Atom> geometry = read_xyz(xyz_name);
	map<string,std::vector<Shell>> basis_set=readBasis(basis_name);
	for (Atom atm: geometry)
	{
		std::vector<Shell> shells = basis_set[atm.element];
		for(const auto&shell:shells)
		{
			auto new_shell = shell;
			new_shell.update_coords(atm.coords);
			basis.push_back(new_shell);
		}
		//basis.insert(basis.end(),shells.begin(),shells.end());
	}
	return basis.size();
}

map<string,std::vector<Shell>> BasisSet::readBasis(string basis_name)
{
	map<string, int> shell2angmom = {{"S", 0}, {"P", 1}, {"D", 2},{"F",3}};
    map<string,std::vector<Shell>> basis_set;
    std::ifstream is(basis_name);
    if (is.good())
    {
      std::string line, rest;
      while (std::getline(is, line) && line != "****") continue;
      bool nextelement = true, nextshell = false;
      std::string cur_element;
      std::vector<Shell> shells;
      // read lines till end
      while (std::getline(is, line))
      {
        // skipping empties and starting with '!' (the comment delimiter)
        if (line.empty() || line[0] == '!') continue;
        if (line == "****")
        {
          if(!shells.empty())
          {
        	  basis_set[cur_element]=shells;
          	  shells.clear();
          }
          nextelement = true;
          nextshell = false;
          continue;
        }
        if (nextelement)
        {
          nextelement = false;
          std::istringstream iss(line);
          iss >> cur_element >> rest;
          transform(cur_element.begin(),cur_element.end(),cur_element.begin(),::toupper);
          nextshell = true;
          continue;
        }
        if (nextshell)
        {
          std::istringstream iss(line);
          std::string shell_label;
          std::size_t n_prims;
          iss >> shell_label >> n_prims >> rest;
            vector<double> exps; vector<double> coeffs;
            for (size_t i=0; i<n_prims; i++)
            {
              while (std::getline(is, line) && (line.empty() || line[0] == '!')) continue;
              std::istringstream iss(line);
              double exponent, coefficient;
              iss >> exponent >> coefficient;
              exps.push_back(exponent);
              coeffs.push_back(coefficient);
            }
            size_t angmom = shell2angmom[shell_label];
            shells.push_back(Shell(angmom,true,exps,coeffs));
        }
      }
}
    return basis_set;
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
	    atoms.push_back(Atom(element_symbol,x,y,z));
	}
	return atoms;


}
