 /*! \file keywords.h
     \brief List of valid keywords
 */
#pragma once
#include <string>
#include <vector>
#include "opencap_exception.h"
#include <map>
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
bool check_keyword(std::string my_key,std::string my_section);
std::map <std::string,std::string> get_params_for_field(std::map<std::string,std::string> params,
		std::string field);


#endif /* INCLUDE_KEYWORDS_H_ */
