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

/*! \file keywords.h
     \brief List of valid keywords
 */
#pragma once

#include <map>
#include <string>
#include <vector>

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
    /** Name of keyword
     */
    std::string name;
    /** Input section keyword belongs to
     */
    std::string input_section;
    /** Valid choices for keyword
     */
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
