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

/*
 * BasisSetParser.cpp
 */

#include "BasisSetParser.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>

#include "opencap_exception.h"
#include "Shell.h"
#include "utils.h"

BasisSetParser::BasisSetParser(std::map<std::string, std::string> parameters)
{

	cart_bf = parameters["cart_bf"];
    transform(cart_bf.begin(),cart_bf.end(),cart_bf.begin(),::toupper);
    for (int i = 0; i < cart_bf.length(); i++)
    	shell2angmom(std::string(1,cart_bf[i]));
	basis_file = parameters["basis_file"];
}

map<string,std::vector<Shell>> BasisSetParser::read_basis()
{
	return read_basis_file();
}

map<string,std::vector<Shell>> BasisSetParser::read_basis_file()
{
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
				  fortran_dfloats_to_efloats(line);
				  std::istringstream iss(line);
				  double exponent, coefficient;
				  iss >> exponent >> coefficient;
				  exps.push_back(exponent);
				  coeffs.push_back(coefficient);
				}
				size_t angmom = shell2angmom(shell_label);
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
				  fortran_dfloats_to_efloats(line);
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
