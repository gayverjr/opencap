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
#include "opencap_exception.h"
#include "keywords.h"
#include "System.h"
#include "utils.h"

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
			if(compare_strings(line,"$geometry"))
				break;
		}
		while (!is.eof())
		{
			std::getline(is, line);
			std::string copied_line = line;
			std::transform(copied_line.begin(), copied_line.end(), copied_line.begin(), ::tolower);
			if (compare_strings(copied_line,"$end"))
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


void parse_section(std::string input_file,std::map<std::string,std::string> &parameters,
		std::string section_name)
{
	std::ifstream is(input_file);
	if (is.good())
	{
	    std::string line, rest;
		while (std::getline(is, line))
		{
			if(compare_strings(line,"$"+section_name))
				break;
		}
		while (!is.eof())
		{
			std::getline(is, line);
			if (compare_strings(line,"$end"))
				break;
			else
			{
				if (line[0] != '!')
				{
					std::istringstream iss(line);
					std::string field,param;
					iss >> field >> param;
					std::transform(field.begin(), field.end(), field.begin(), ::tolower);
					if (check_keyword(field,section_name,param))
						parameters[field] = param;
					else
						opencap_throw("Invalid keyword in " + section_name +" section:" + field );
				}
			}
			}
		}
}

System get_System(std::string input_file, std::map<std::string,std::string> params)
{
	//first check that we've got what we need to at least try to construct a system
	if (params.find("molecule")==params.end())
		opencap_throw("Missing required keyword: molecule. Please choose one of the following: qchem_fchk, "
				"molcas_rassi,read, molden");
	if(params.find("basis_file")==params.end())
		opencap_throw("Error: Need to specify a basis set file using the basis_file keyword.");
	//geometry
	if(compare_strings(params["molecule"],"inline"))
		return System(parse_geometry(input_file),params);
	else
		return System(params["basis_file"],params["molecule"]);

}

std::tuple<System,std::map<std::string,std::string>> parse_input(std::string input_file)
{
	//TODO: check that fields are present
	std::map<std::string,std::string> parameters;
	System my_sys;
	try
	{
		//parse job section
		parse_section(input_file,parameters,"job");
		//parse system
		parse_section(input_file,parameters,"system");
		if(compare_strings(parameters["jobtype"],"projected_cap"))
		{
			//parse cap_parameters
			parse_section(input_file,parameters,"projected_cap");
		}
		else
			opencap_throw("Invalid jobtype: " + parameters["jobtype"]);
		my_sys = get_System(input_file,get_params_for_field(parameters,"system"));
	}
	catch (exception &e)
	{
		opencap_rethrow("Failed to parse input file.");
	}
	return std::make_tuple(my_sys,parameters);
}


