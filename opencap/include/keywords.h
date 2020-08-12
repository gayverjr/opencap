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

/*! \brief Struct denoting a keyword which can have many different values
 *
 *
 */
struct open_key
{
	/** Name of keyword
	  */
	std::string name;
	/** Input section keyword belongs to
	  */
	std::string input_section;
};

/*! \brief Struct denoting a keyword which can only have certain values
 *
 *
 */
struct closed_key
{
	std::string name;
	std::string input_section;
	std::vector<std::string> valid_options;
};

/** Checks whether a keyword is valid.
*/
bool check_keyword(std::string my_key,std::string my_section,std::string my_value);
/** Collects parameters which apply to input field.
 */
std::map <std::string,std::string> get_params_for_field(std::map<std::string,std::string> params,
		std::string field);


#endif /* INCLUDE_KEYWORDS_H_ */
