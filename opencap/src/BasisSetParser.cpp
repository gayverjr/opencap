/*
 * BasisSetParser.cpp
 *
 *  Created on: Apr 6, 2020
 *      Author: JG
 */
#include "BasisSetParser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <map>
#include <list>
#include "Shell.h"

BasisSetParser::BasisSetParser(std::map<std::string, std::string> parameters)
{
	cart_bf = parameters["cart_bf"];
	//change cart_bf to upper case
    transform(cart_bf.begin(),cart_bf.end(),cart_bf.begin(),::toupper);
	basis_file = parameters["basis_file"];
}

map<string,std::vector<Shell>> BasisSetParser::read_basis()
{
	return read_basis_file();
}

map<string,std::vector<Shell>> BasisSetParser::read_basis_file()
{
	map<string, int> shell2angmom = {{"S", 0}, {"P", 1}, {"D", 2},{"F",3},{"G",4},{"H",5}};
    map<string,std::vector<Shell>> basis_set;
    std::ifstream is(basis_file);
    if (is.good())
    {
      std::string line, rest;
      while (std::getline(is, line) && line != "****") continue;
      bool nextelement = true, nextshell = false;
      std::string cur_element;
      std::vector<Shell> shells;
      // read lines till end
      while (std::getline(is, line))
      {
        // skipping empties and starting with '!' (the comment delimiter)
        if (line.empty() || line[0] == '!') continue;
        if (line == "****")
        {
          if(!shells.empty())
          {
        	  basis_set[cur_element]=shells;
          	  shells.clear();
          }
          nextelement = true;
          nextshell = false;
          continue;
        }
        if (nextelement)
        {
          nextelement = false;
          std::istringstream iss(line);
          iss >> cur_element >> rest;
          transform(cur_element.begin(),cur_element.end(),cur_element.begin(),::tolower);
          cur_element[0] = toupper(cur_element[0]);
          nextshell = true;
          continue;
        }
        if (nextshell)
        {

          std::istringstream iss(line);
          std::string shell_label;
          std::size_t n_prims;
          iss >> shell_label >> n_prims >> rest;
          if (shell_label!="SP")
          {
				vector<double> exps,coeffs;
				for (size_t i=0; i<n_prims; i++)
				{
				  while (std::getline(is, line) && (line.empty() || line[0] == '!')) continue;
				  std::istringstream iss(line);
				  double exponent, coefficient;
				  iss >> exponent >> coefficient;
				  exps.push_back(exponent);
				  coeffs.push_back(coefficient);
				}
				size_t angmom = shell2angmom[shell_label];
				if (cart_bf.find(shell_label) != std::string::npos)
					shells.push_back(Shell(angmom,false,exps,coeffs));
				else
					shells.push_back(Shell(angmom,true,exps,coeffs));
          }
          if (shell_label=="SP")
          {
				vector<double> exps, s_coeffs, p_coeffs;
				for (size_t i=0; i<n_prims; i++)
				{
				  while (std::getline(is, line) && (line.empty() || line[0] == '!')) continue;
				  std::istringstream iss(line);
				  double exponent, s_coeff, p_coeff;
				  iss >> exponent >> s_coeff >> p_coeff;
				  exps.push_back(exponent);
				  s_coeffs.push_back(s_coeff);
				  p_coeffs.push_back(p_coeff);
				}
				if (cart_bf.find("S") != std::string::npos)
					shells.push_back(Shell(0,false,exps,s_coeffs));
				else
					shells.push_back(Shell(0,true,exps,s_coeffs));
				if (cart_bf.find("P") != std::string::npos)
					shells.push_back(Shell(1,false,exps,p_coeffs));
				else
					shells.push_back(Shell(1,true,exps,p_coeffs));
          }
        }
      }
}
    return basis_set;
}
