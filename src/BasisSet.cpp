#include <vector>
#include "BasisSet.h"
#include "Atom.h"
#include <math.h>
#include <algorithm>
#include <string>
#include <map>
#include <list>
#include "Shell.h"
#include <array>
#include <limits>
#include <iostream>
#include "BasisSetParser.h"


BasisSet::BasisSet(std::vector<Atom> geometry,std::map<std::string, std::string> parameters)
{
	//TODO don't assume there's a basis name
	std::string basis_name = parameters["basis_file"];
	std::string cart_bf = parameters["cart_bf"];
	BasisSetParser parser (parameters);
	build_basis_set(geometry,parser.read_basis());
	Nshells = basis.size();
	Nbasis = calc_basis_size();
	std::cout << "Number of basis functions:" << Nbasis << std::endl;
	name = "gen";
}

BasisSet::BasisSet()
{
	Nshells=0;
	Nbasis=0;
	name="undefined";
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

void BasisSet::build_basis_set(std::vector<Atom> geometry,map<string,std::vector<Shell>> all_shells)
{

	for (Atom atm: geometry)
	{
		std::vector<Shell> my_shells = all_shells[atm.symbol];
		for(const auto&shell:my_shells)
		{
			auto new_shell = shell;
			new_shell.update_coords(atm.coords);
			basis.push_back(new_shell);
		}
	}
	Nshells = basis.size();
	Nbasis = calc_basis_size();
	name = "user defined";
}

int BasisSet::max_L()
{
	int l_max = 0;
	for(const auto&shell:basis)
	{
		if(shell.l > l_max)
			l_max = shell.l;
	}
	return l_max;
}

double BasisSet::alpha_max(Atom atm)
{
	double max_val = 0;
	std::vector<Shell> shells = shells_on_center(atm);
	for(const auto&shell:shells)
	{
		for(const double exp: shell.exps)
		{
			if(exp>max_val)
				max_val=exp;
		}
	}
	return max_val;
}

std::vector<double> BasisSet::alpha_min(Atom atm)
{
	std::vector<Shell> shells = shells_on_center(atm);
	std::vector<double> min_alpha(max_L()+1,0.0);
	for(const auto&shell:shells)
	{
		for(const double exp: shell.exps)
		{
			if(exp<min_alpha[shell.l] || min_alpha[shell.l]==0.0)
				min_alpha[shell.l]=exp;
		}
	}
	return min_alpha;
}

std::vector<Shell> BasisSet::shells_on_center(Atom atm)
{
	std::vector<Shell> shells;
	for(const auto&shell:basis)
	{
		if(shell.origin[0]==atm.coords[0] && shell.origin[1]==atm.coords[1] && shell.origin[2]==atm.coords[2])
			shells.push_back(shell);
	}
	return shells;
}


