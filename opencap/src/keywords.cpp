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
#include <algorithm>
#include <string>
#include <map>
#include <iostream>

bool check_keyword(std::string my_key,std::string my_section)
{
	std::transform(my_key.begin(), my_key.end(), my_key.begin(), ::tolower);
	std::vector<valid_key> valid_keys
	{
		{"title","job"},
		{"jobtype","job"},
		{"bohr_coordinates","system"},
		{"cart_bf","system"},
		{"basis_file","system"},
		{"molecule","system"},
		{"method","projected_cap"},
		{"package","projected_cap"},
		{"nstates","projected_cap"},
		{"molcas_output","projected_cap"},
		{"rassi_h5","projected_cap"},
		{"h0_file","projected_cap"},
		{"qchem_output","projected_cap"},
		{"qchem_fchk","projected_cap"},
		{"cap_type","projected_cap"},
		{"cap_x","projected_cap"},
		{"cap_y","projected_cap"},
		{"cap_z","projected_cap"},
		{"r_cut","projected_cap"},
		{"radial_precision","projected_cap"},
		{"angular_points","projected_cap"}
	};
	for(auto key:valid_keys)
	{
		if(key.name==my_key && key.input_section==my_section)
			return true;
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
		if(check_keyword(pair.first,field))
			params_for_field[pair.first]=pair.second;
	}
	return params_for_field;
}



