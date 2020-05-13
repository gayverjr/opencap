/*
 * keywords.h
 *
 *  Created on: May 13, 2020
 *      Author: JG
 */
#include <string>
#include <vector>
#include "opencap_exception.h"
#ifndef INCLUDE_KEYWORDS_H_
#define INCLUDE_KEYWORDS_H_

struct valid_key
{
	std::string name;
	std::string input_section;
};

bool check_keyword(std::string my_key,std::string my_section)
{
	std::vector<valid_key> valid_keys
	{
		{"title","job"},
		{"jobtype","job"},
		{"bohr_coordinates","system"},
		{"cart_bf","system"},
		{"method","system"},
		{"package","system"},
		{"rassi_h5","system"},
		{"nstates","system"},
		{"molcas_output","system"},
		{"rassi_h5","system"},
		{"h0_file","system"},
		{"basis_file","system"},
		{"qc_output","system"},
		{"fchk_file","system"},
		{"cap_type","cap_parameters"},
		{"cap_x","cap_parameters"},
		{"cap_y","cap_parameters"},
		{"cap_z","cap_parameters"},
		{"r_cut","cap_parameters"},
		{"radial_precision","cap_parameters"},
		{"angular_points","cap_parameters"}
	};
	for(auto key:valid_keys)
	{
		if(key.name==my_key && key.input_section==my_section)
			return true;
	}
	return false;
};


#endif /* INCLUDE_KEYWORDS_H_ */
