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
 * keywords.cpp
 */

#include "keywords.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "opencap_exception.h"
#include "utils.h"

bool check_keyword(std::string my_key,std::string my_section,std::string my_value)
{
	std::transform(my_key.begin(), my_key.end(), my_key.begin(), ::tolower);
	std::vector<open_key> open_keys
	{
		{"eta_step","trajectory"},
		{"nsteps","trajectory"},
		{"save_trajectory","trajectory"},
		{"basis_file","system"},
		{"method","projected_cap"},
		{"nstates","projected_cap"},
		{"molcas_output","projected_cap"},
		{"rassi_h5","projected_cap"},
		{"h0_file","projected_cap"},
		{"qchem_output","projected_cap"},
		{"qchem_fchk","projected_cap"},
		{"cap_x","projected_cap"},
		{"cap_y","projected_cap"},
		{"cap_z","projected_cap"},
		{"r_cut","projected_cap"},
		{"radial_precision","projected_cap"},
		{"angular_points","projected_cap"},
		{"cart_bf","system"}
	};
	std::vector<closed_key> closed_keys
	{
		{"package","projected_cap",{"openmolcas","qchem"}},
		{"cap_type","projected_cap",{"box","voronoi"}},
		{"molecule","system",{"inline","molden","qchem_fchk","molcas_rassi"}},
		{"bohr_coordinates","system",{"true","false"}},
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



