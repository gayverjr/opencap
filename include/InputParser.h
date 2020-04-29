/*
 * InputParser.h
 *
 *  Created on: Apr 1, 2020
 *      Author: JG
 */

#ifndef INPUTPARSER_H_
#define INPUTPARSER_H_
#include <map>
#include <tuple>
#include <string>


std::tuple<std::vector<Atom>,std::map<std::string,std::string>> parse_input(std::string input_file);
std::vector<Atom> parse_geometry(std::string input_file);
void parse_system_info(std::string input_file,std::map<std::string,std::string> &parameters);
void parse_cap_parameters(std::string input_file,std::map<std::string,std::string> &parameters);



#endif /* INPUTPARSER_H_ */
