/*
 * molden_parser.cpp
 *
 *  Created on: Jul 14, 2020
 *      Author: JG
 */
#include "molden_parser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "utils.h"
#include "Atom.h"
#include "BasisSet.h"
#include "Shell.h"
#include <map>


std::vector<Atom> read_geometry_from_molden(std::string filename)
{
	std::vector<Atom> atoms;
	std::ifstream is(filename);
	if(is.good())
	{
		//first get the number of atoms
    	std::string line, rest;
    	std::getline(is, line);
		while(line.find("N_ATOMS")== std::string::npos)
			std::getline(is,line);
    	std::getline(is, line);
		size_t num_atoms = stoi(split(line,' ').back());
		//now lets read in the atoms
		while(line.find("ATOMS")== std::string::npos)
			std::getline(is,line);
		for(size_t i=1;i<=num_atoms;i++)
		{
			std::getline(is,line);
			size_t charge  = stoi(split(line,' ')[2]);
			double x_coord = stod(split(line,' ')[3]);
			double y_coord = stod(split(line,' ')[4]);
			double z_coord = stod(split(line,' ')[5]);
			atoms.push_back(Atom(charge,x_coord,y_coord,z_coord));
		}
	}
	return atoms;
}

Shell read_shell_from_molden(std::string line,std::ifstream &is,std::array<double,3> cur_coords)
{
	map<string, int> shell2angmom = {{"S", 0}, {"P", 1}, {"D", 2},{"F",3},{"G",4},{"H",5}};
    std::istringstream iss(line);
    std::string rest;
    std::string shell_label;
    std::size_t n_prims;
    iss >> shell_label >> n_prims >> rest;
    //make shell label upper case
    transform(shell_label.begin(),shell_label.end(),shell_label.begin(),::toupper);
    size_t angmom = shell2angmom[shell_label];
    Shell my_shell(angmom,cur_coords);
    for(size_t i=1;i<=n_prims;i++)
    {
		std::getline(is,line);
	    std::istringstream iss_prim(line);
		double exp,coeff;
		iss_prim >> exp >> coeff;
		my_shell.add_primitive(exp,coeff);
    }
    return my_shell;
}

BasisSet read_basis_from_molden(std::string filename,std::vector<Atom> atoms)
{
	std::ifstream is(filename);
	bool harmonic_d = false;
	bool harmonic_f = false;
	bool harmonic_g = false;
	BasisSet bs;
	for(auto atm:atoms)
		bs.centers.push_back(atm.coords);
	if(is.good())
	{
    	std::string line, rest;
		//first let's set flags for 5d,7f,9G
		while(is.peek() != EOF )
		{
	    	std::getline(is, line);
	    	if(line.find("5D")!= std::string::npos)
	    		harmonic_d = true;
	    	if(line.find("7F")!= std::string::npos)
	    		harmonic_f = true;
	    	if(line.find("9G")!= std::string::npos)
	    		harmonic_g = true;
		}
		//back to the beginning
		is.seekg (0, ios::beg);
    	std::getline(is, line);
		while(line.find("GTO")== std::string::npos)
			std::getline(is,line);
		//read til we hit MO section
		std::getline(is,line);
		while(line.find("MO")==std::string::npos)
		{
			size_t atm_idx;
			//case 1: atom index
			if(split(line,' ').size()==1)
			    atm_idx = stoi(split(line,' ')[0])-1;
			//case 2: shell and number of primitives
			else if(split(line,' ').size()==2)
			{
				Shell my_shell = read_shell_from_molden(line,is,atoms[atm_idx].coords);
				//now lets check for pure/cartesian flags
				if(my_shell.l==2 && !harmonic_d)
					my_shell.pure=false;
				else if(my_shell.l==3 && !harmonic_f)
					my_shell.pure=false;
				else if(my_shell.l==4 && !harmonic_g)
					my_shell.pure=false;
				bs.add_shell(my_shell);
			}
			std::getline(is,line);
		}
	}
	bs.normalize();
	return bs;
}
