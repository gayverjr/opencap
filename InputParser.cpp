/*
 * InputParser.cpp
 *
 *  Created on: Apr 1, 2020
 *      Author: JG
 */

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <map>
#include <list>
#include <fstream>
#include <tuple>
#include <locale>
#include "Atom.h"
#include "InputParser.h"

void parse_system_info(std::string input_file,std::map<std::string,std::string> &parameters)
{
	std::vector<std::string> valid_fields = {"bohr_coordinates","cart_bf","package","method",
											"nstates","basis_set","basis_file",
											"fchk_file","qc_output","rassi_h5","molcas_output"};
	std::vector<Atom> atoms;
	std::ifstream is(input_file);
	if (is.good())
	{
	    std::string line, rest;
		while (std::getline(is, line))
		{
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			if(line=="$system")
				break;
		}
		while (!is.eof())
		{
			std::getline(is, line);
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			if (line=="$end")
				break;
			else
			{
				if (line[0] != '!')
				{
					std::istringstream iss(line);
					std::string field,param;
					iss >> field >> param;
					std::transform(field.begin(), field.end(), field.begin(), ::tolower);
					if (std::find(valid_fields.begin(), valid_fields.end(), field) != valid_fields.end())
					{
						parameters[field] = param;
					}
					else
						std::cout << "Invalid field:" << field << std::endl;
					}
				}
			}
		}
}

void parse_cap_params(std::string input_file,std::map<std::string,std::string> &parameters)
{
	std::vector<std::string> valid_fields = {"cap_type","cap_x","cap_y","cap_z","r_cut",
			"radial_precision","angular_points"};
	std::vector<Atom> atoms;
	std::ifstream is(input_file);
	if (is.good())
	{
	    std::string line, rest;
		while (std::getline(is, line))
		{
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			if(line=="$cap_parameters")
				break;
		}
		while (!is.eof())
		{
			std::getline(is, line);
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			if (line=="$end")
				break;
			else
			{
				if (line[0] != '!')
				{
					std::istringstream iss(line);
					std::string field,param;
					iss >> field >> param;
					std::transform(field.begin(), field.end(), field.begin(), ::tolower);
					if (std::find(valid_fields.begin(), valid_fields.end(), field) != valid_fields.end())
					{
						parameters[field] = param;
					}
					else
						std::cout << "Invalid field:" << field << std::endl;
					}
				}
			}
		}
}

std::vector<Atom> parse_geometry(std::string input_file)
{
	std::vector<Atom> atoms;
	std::ifstream is(input_file);
	if (is.good())
	{
	    std::string line, rest;
		while (std::getline(is, line))
		{
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			if(line=="$geometry")
				break;
		}
		while (!is.eof())
		{
			std::getline(is, line);
			std::string copied_line = line;
			std::transform(copied_line.begin(), copied_line.end(), copied_line.begin(), ::tolower);
			if (copied_line=="$end")
				break;
			else
			{
					std::istringstream iss(line);
					std::string element_symbol;
					double x, y, z;
					iss >> element_symbol >> x >> y >> z;
					transform(element_symbol.begin(),element_symbol.end(),element_symbol.begin(),::tolower);
					atoms.push_back(Atom(element_symbol,x,y,z));
			}
		}
	}
	return atoms;
}

std::tuple<std::vector<Atom>,std::map<std::string,std::string>> parse_input(std::string input_file)
{
	std::map<std::string,std::string> parameters;
	//run parser for each input field
	std::vector<Atom>atoms = parse_geometry(input_file);
	parse_cap_params(input_file,parameters);
	parse_system_info(input_file,parameters);
	auto inp_data = std::make_tuple(atoms, parameters);
	return inp_data;
}


