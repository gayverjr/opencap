/*
 * keywords.cpp
 *
 *  Created on: May 29, 2020
 *      Author: JG
 */
#include <string>
#include <vector>
#include "opencap_exception.h"
#include "keywords.h"
#include "utils.h"
#include <algorithm>
#include <string>
#include <map>
#include <iostream>

bool check_keyword(std::string my_key,std::string my_section,std::string my_value)
{
	std::transform(my_key.begin(), my_key.end(), my_key.begin(), ::tolower);
	std::vector<open_key> open_keys
	{
		{"title","job"},
		{"jobtype","job"},
		{"basis_file","system"},
		{"method","perturb_cap"},
		{"nstates","perturb_cap"},
		{"molcas_output","perturb_cap"},
		{"rassi_h5","perturb_cap"},
		{"h0_file","perturb_cap"},
		{"qchem_output","perturb_cap"},
		{"qchem_fchk","perturb_cap"},
		{"cap_x","perturb_cap"},
		{"cap_y","perturb_cap"},
		{"cap_z","perturb_cap"},
		{"r_cut","perturb_cap"},
		{"radial_precision","perturb_cap"},
		{"angular_points","perturb_cap"},
	};
	std::vector<closed_key> closed_keys
	{
		{"package","perturb_cap",{"openmolcas","qchem"}},
		{"cap_type","perturb_cap",{"box","voronoi"}},
		{"molecule","system",{"inline","molden","qchem_fchk","molcas_rassi"}},
		{"bohr_coordinates","system",{"true","false"}},
		{"cart_bf","system",{"d","df","dfg","f","fg","g"}},
	};

	for(auto key:open_keys)
	{
		if(compare_strings(key.name,my_key) && compare_strings(key.input_section,my_section))
			return true;
	}
	for(auto key:closed_keys)
	{
		if(compare_strings(key.name,my_key))
		{
			std::string val_copy = my_value;
			std::transform(val_copy.begin(), val_copy.end(), val_copy.begin(), ::tolower);
			if(std::find(key.valid_options.begin(), key.valid_options.end(), val_copy) != key.valid_options.end())
				return true;
			else
				opencap_throw("Invalid choice: \'"+ my_value + "\' for keyword: \'" + my_key + "\'");
		}
	}
	return false;
}

std::map <std::string,std::string> get_params_for_field(
		std::map<std::string,std::string> params,
		std::string field)
{
	std::map<std::string,std::string> params_for_field;
	for (const auto &pair : params)
	{
		if(check_keyword(pair.first,field,pair.second))
			params_for_field[pair.first]=pair.second;
	}
	return params_for_field;
}



