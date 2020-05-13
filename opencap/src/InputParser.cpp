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


void parse_section(std::string input_file,std::map<std::string,std::string> &parameters,
		std::string section_name)
{
	std::ifstream is(input_file);
	if (is.good())
	{
	    std::string line, rest;
		while (std::getline(is, line))
		{
			std::transform(line.begin(), line.end(), line.begin(), ::tolower);
			if(line=="$"+section_name)
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
					if (check_keyword(field,section_name))
						parameters[field] = param;
					else
						opencap_throw("Invalid keyword in " + section_name +" section:" + field );
				}
			}
			}
		}
}

std::tuple<std::vector<Atom>,std::map<std::string,std::string>> parse_input(std::string input_file)
{
	//TODO: check that fields are present
	std::map<std::string,std::string> parameters;
	std::vector<Atom>atoms;
	try
	{
		//run parser for each input field
		atoms = parse_geometry(input_file);
		//parse job section
		parse_section(input_file,parameters,"job");
		if(parameters["jobtype"]=="projected_cap")
		{
			//parse system
			parse_section(input_file,parameters,"system");
			//parse cap_parameters
			parse_section(input_file,parameters,"projected_cap");
		}
		else
			opencap_throw("Invalid jobtype: " + parameters["jobtype"]);

	}
	catch (exception &e)
	{
		opencap_rethrow("Failed to parse input file.");
	}
	return std::make_tuple(atoms,parameters);
}


