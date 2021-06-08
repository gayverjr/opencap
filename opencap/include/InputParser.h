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

/*! \file InputParser.h
     \brief Class for parsing OpenCAP input files.
 */
#ifndef INPUTPARSER_H_
#define INPUTPARSER_H_

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "Atom.h"
#include "System.h"



/** Parses an input file.
 * \param input_file: OpenCAP input file
 * \return A tuple, 1st element is the geometry, second element is the parameters map
*/
std::tuple<System,std::map<std::string,std::string>> parse_input(std::string input_file);
/** Parses the geometry section.
 * \param input_file: OpenCAP input file
 * \return Geometry of molecular system
 */
std::vector<Atom> parse_geometry(std::string input_file);
/** Parses a section of the input file.
 * \param input_file: OpenCAP input file
 * \param parameters: Parameters map
 * \param section_name: Name of section
 */
void parse_section(std::string input_file,std::map<std::string,std::string> &parameters,
		std::string section_name);

#endif /* INPUTPARSER_H_ */
