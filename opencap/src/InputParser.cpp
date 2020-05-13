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

void parse_system_info(std::string input_file,std::map<std::string,std::string> &parameters)
{
	std::vector<std::string> valid_fields = {"bohr_coordinates","cart_bf","package","method",
											"nstates","basis_set","basis_file","h0_file",
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
	check_parameters(parameters);
	auto inp_data = std::make_tuple(atoms, parameters);
	return inp_data;
}

void verify_cap_parameters(std::map<std::string,std::string> parameters)
{
	std::vector<std::string> missing_keys;
	if(parameters.find("cap_type")==parameters.end())
		opencap_throw("Error: Missing cap_type keyword.");
	if(parameters["cap_type"]=="box")
	{
		if(parameters.find("cap_x")==parameters.end())
			missing_keys.push_back("cap_x");
		if(parameters.find("cap_y")==parameters.end())
			missing_keys.push_back("cap_y");
		if (parameters.find("cap_z")==parameters.end())
			missing_keys.push_back("cap_z");
	}
	else if (parameters["cap_type"]=="voronoi")
	{
		if(parameters.find("r_cut")==parameters.end())
			missing_keys.push_back("r_cut");
	}
	else
		opencap_throw("Error: only box and voronoi CAPs supported.");
	if(missing_keys.size()!=0)
	{
		std::string error_str = "Missing CAP keywords:";
		for (auto key: missing_keys)
			error_str+=key+",";
		opencap_throw(error_str);
	}
}

void verify_method(std::map<std::string,std::string> parameters,std::string key)
{
	std::string package_name = parameters[key];
	if(parameters.find("method")==parameters.end())
	{
		std::cout << "Error: missing the 'method' keyword. Please choose a supported package/method." << std::endl;
	}
	std::string method = parameters["method"];
	if (package_name=="qchem")
	{
		std::vector<std::string> supported = {"eomea","eomee","eomip"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
		{
			std::cout << "Error: unsupported Q-Chem method. OpenCAP currently supports: 'eom-ea-ccsd','eom-ee-ccsd','eom-ip-ccsd'."<< std::endl;
		}
		if (parameters.find("fchk_file")==parameters.end())
		{
			std::cout << "Error: missing keyword: fchk_file." << std::endl;
		}
		if (parameters.find("qc_output")==parameters.end() && parameters.find("h0_file")==parameters.end())
		{
			std::cout << "Error: Need to specify H0 Hamiltonian via \"qc_output\" or \"h0_file\" fields." << std::endl;
		}
	}
	else if(package_name=="openmolcas")
	{
		std::vector<std::string> supported = {"ms-caspt2","xms-caspt2","pc-nevpt2","sc-nevpt2"};
		if (std::find(supported.begin(), supported.end(), method) == supported.end())
		{
			std::cout << "Error: unsupported OpenMolcas method. OpenCAP currently supports: 'dmrgscf','nevpt2','caspt2','rasscf'."<< std::endl;
		}
		if (parameters.find("rassi_h5")==parameters.end())
		{
			std::cout << "Error: missing keyword: rassi_h5." << std::endl;
		}
		if (parameters.find("molcas_output")==parameters.end() && parameters.find("h0_file")==parameters.end())
		{
			std::cout << "Error: Need to specify H0 Hamiltonian via \"molcas_output\" or \"h0_file\" fields." << std::endl;
		}
	}
	else
	{
		std::cout << "Error: unsupported package. Only Q-Chem and OpenMolcas are currently supported." << std::endl;
	}
}

//check that requirements have been specified, and set optional params to their defaults
void check_parameters(std::map<std::string,std::string> parameters)
{
	std::map<std::string, std::string> defaults =
	{{"bohr_coordinates", "false"}, {"radial_precision", "14"}, {"angular_points", "590"},{"cart_bf",""},{"basis_file",""}};
	std::vector<std::string> requirements = {"package","nstates","cap_type","basis_set"};
	//first lets verify that all of our requirements are there
	for (std::string key:requirements)
	{
		if (parameters.find(key)!=parameters.end())
		{
			if (key=="cap_type")
			{
				try
				{
					verify_cap_parameters(parameters);
				}
				catch (exception &e)
				{
					opencap_rethrow("$cap_parameters section is invalid.");
				}
			}
			if (key=="basis_set" && parameters[key]=="gen")
			{
				if (parameters.find("basis_file")==parameters.end())
				{
					std::cout << "Need to specify a basis file when using gen basis." << std::endl;
				}
			}
			else if(key=="package")
			{
				if (parameters.find("h0_file")==parameters.end())
				{
					verify_method(parameters,key);
				}
			}
		}
		else
		{
			std::cout << "Error: missing required key \"" << key << "\""<< std::endl;
		}
	}
	//now let's set optional parameters to their defaults if not specified
	for (const auto &pair:defaults)
	{
		if(parameters.find(pair.first)==parameters.end())
			parameters[pair.first]=pair.second;
	}
}

