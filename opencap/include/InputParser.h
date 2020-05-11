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
/** Parses the system section.
 */
void parse_system_info(std::string input_file,std::map<std::string,std::string> &parameters);
/** Parses the %CAP section.
 */
void parse_cap_parameters(std::string input_file,std::map<std::string,std::string> &parameters);



#endif /* INPUTPARSER_H_ */
