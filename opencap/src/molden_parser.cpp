/*Copyright (c) 2021 James Gayvert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
 * molden_parser.cpp
 *
 */

#include "molden_parser.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "Atom.h"
#include "BasisSet.h"
#include "opencap_exception.h"
#include "Shell.h"
#include "utils.h"


std::vector<Atom> read_geometry_from_molden(std::string filename)
{
	std::vector<Atom> atoms;
	size_t num_atoms = 0;
	bool angs = false;
	bool n_atoms_found = true;
	std::ifstream is(filename);
	if(is.good())
	{
		//first get the number of atoms
    	std::string line, rest;
    	std::getline(is, line);
		while(line.find("N_ATOMS")== std::string::npos && (line.find("N_Atoms")== std::string::npos))
		{
    		if (is.peek()==EOF)
    		{
    			n_atoms_found = false;
    			//back to beginning
    			is.seekg (0, ios::beg);
    			break;
    		}
    		else
    			std::getline(is, line);
		}
	   	std::getline(is, line);
		if(n_atoms_found)
			num_atoms = stoi(split(line,' ').back());
		//now lets read in the atoms
		while(line.find("ATOMS")== std::string::npos && line.find("Atoms")== std::string::npos)
		{
			std::getline(is,line);
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before \"ATOMS\" section was found.");
		}
		if(line.find("Angs")!=std::string::npos || line.find("ANGS")!=std::string::npos)
			angs = true;
		std::getline(is,line);
		while(line.find("[")==std::string::npos)
		{
    		if (is.peek()==EOF)
    			opencap_throw("Error: Reached end of file before the atoms could be parsed.");
			size_t charge  = stoi(split(line,' ')[2]);
			double x_coord = stod(split(line,' ')[3]);
			double y_coord = stod(split(line,' ')[4]);
			double z_coord = stod(split(line,' ')[5]);
			Atom atm(charge,x_coord,y_coord,z_coord);
			if(angs)
				atm.ang_to_bohr();
			atoms.push_back(atm);
			std::getline(is,line);
		}
	}
	if(n_atoms_found && num_atoms!=atoms.size())
	{
		opencap_throw("Error: Number of atoms found does not match those specified in \"N_ATOMS\" section.");
	}
	else if(atoms.size()==0)
	{
		opencap_throw("Error: No atoms found.");
	}
	return atoms;
}

Shell read_shell_from_molden(std::string line,std::ifstream &is,std::array<double,3> cur_coords)
{
    std::istringstream iss(line);
    std::string rest;
    std::string shell_label;
    std::size_t n_prims;
    iss >> shell_label >> n_prims >> rest;
    //make shell label upper case
    transform(shell_label.begin(),shell_label.end(),shell_label.begin(),::toupper);
    size_t angmom = shell2angmom(shell_label);
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
	    	if(line.find("5D")!= std::string::npos || line.find("5d")!= std::string::npos)
	    		harmonic_d = true;
	    	if(line.find("7F")!= std::string::npos || line.find("7f")!= std::string::npos)
	    		harmonic_f = true;
	    	if(line.find("9G")!= std::string::npos || line.find("9g")!= std::string::npos)
	    		harmonic_g = true;
		}
		//back to the beginning
		is.seekg (0, ios::beg);
    	std::getline(is, line);
		while(line.find("GTO")== std::string::npos)
			std::getline(is,line);
		std::getline(is,line);
		//read til we hit the next section
		while(line.find("[")==std::string::npos)
		{
			size_t atm_idx;
			//case 1: atom index
			if(split(line,' ').size()>0 && is_number(split(line,' ')[0]))
			    atm_idx = stoi(split(line,' ')[0])-1;
			//case 2: shell and number of primitives
			else if(split(line,' ').size()>0 && is_letter(split(line,' ')[0]))
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
