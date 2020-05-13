 /*! \file keywords.h
     \brief List of valid keywords
 */
#include <string>
#include <vector>
#include "opencap_exception.h"
#ifndef INCLUDE_KEYWORDS_H_
#define INCLUDE_KEYWORDS_H_

/*! \brief Struct denoting a valid keyword
 *
 *
 */
struct valid_key
{
	/** Name of keyword
	  */
	std::string name;
	/** Input section keyword belongs to
	  */
	std::string input_section;
};

/** Checks whether a keyword is valid.
*/
bool check_keyword(std::string my_key,std::string my_section)
{
	std::vector<valid_key> valid_keys
	{
		{"title","job"},
		{"jobtype","job"},
		{"bohr_coordinates","system"},
		{"cart_bf","system"},
		{"basis_file","system"},
		{"method","projected_cap"},
		{"package","projected_cap"},
		{"rassi_h5","projected_cap"},
		{"nstates","projected_cap"},
		{"molcas_output","projected_cap"},
		{"rassi_h5","projected_cap"},
		{"h0_file","projected_cap"},
		{"qc_output","projected_cap"},
		{"fchk_file","projected_cap"},
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
};


#endif /* INCLUDE_KEYWORDS_H_ */
