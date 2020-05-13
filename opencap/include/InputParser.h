 /*! \file InputParser.h
     \brief Class for parsing OpenCAP input files.
 */
#ifndef INPUTPARSER_H_
#define INPUTPARSER_H_
#include <map>
#include <tuple>
#include <string>


/** Parses an input file.
*/
std::tuple<std::vector<Atom>,std::map<std::string,std::string>> parse_input(std::string input_file);
/** Parses the geometry section.
 */
std::vector<Atom> parse_geometry(std::string input_file);
void parse_section(std::string input_file,std::map<std::string,std::string> &parameters,
		std::string section_name);

#endif /* INPUTPARSER_H_ */
