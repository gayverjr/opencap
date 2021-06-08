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

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "opencap_exception.h"
#include "Shell.h"

static map<string, int> angmom_map = {{"S", 0}, {"P", 1}, {"D", 2},{"F",3},{"G",4}};
double fact2(int n)
{
    int res = 1;
    for (int i=n; i>=0; i=i-2)
    {
        if (i==0 || i==1)
            return res;
        else
            res *= i;
    }
    return static_cast<double>(res);
}

int parity(int i)
{
  return i%2 ? -1 : 1;
}

double fact(int n)
{
    int res = 1;
    for(int i = 1; i <= n; ++i)
        res *= i;
    return static_cast<double>(res);
}

double binom(int p, int q)
{
	if(q<0 || q>p)
		return 0;
	return fact(p)/fact(q)/fact(p-q);
}


std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
	  if (!token.empty())
		  tokens.push_back(token);
   }
   return tokens;
}


void fortran_dfloats_to_efloats(std::string& str)
{
	for(auto& ch: str)
	{
		if (ch == 'd') ch = 'e';
		if (ch == 'D') ch = 'E';
	}
}

int shell2angmom(std::string shell_label)
{
	if(angmom_map.find(shell_label)==angmom_map.end())
		opencap_throw("Error: Only up to G type shells are supported.");
	return angmom_map[shell_label];
}

bool compare_strings(std::string s1, std::string s2)
{
	std::transform(s1.begin(), s1.end(), s1.begin(), ::tolower);
	std::transform(s2.begin(), s2.end(), s2.begin(), ::tolower);
	return s1==s2;
}

bool is_number(const std::string &s)
{
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

bool is_letter(const std::string &s)
{
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isalpha);
}


Eigen::MatrixXd read_matrix(size_t N, std::string filename)
{
	Eigen::MatrixXd h0(N,N);
	h0= Eigen::MatrixXd::Zero(N,N);
	std::ifstream is(filename);
	if (is.good())
	{
		//first line should be diagonal or full
		std::string line;
		std::getline(is,line);
		std::string mat_type = line;
		std::transform(mat_type.begin(), mat_type.end(), mat_type.begin(), ::tolower);
		if (compare_strings(mat_type,"diagonal"))
		{
			for(size_t i=0;i<N;i++)
			{
				std::getline(is,line);
				h0(i,i)=std::stod(line);
			}
		}
		else if (compare_strings(mat_type,"full"))
		{
			for (size_t i=0;i<N;i++)
			{
				std::getline(is,line);
				std::vector<std::string> tokens = split(line,' ');
				for(size_t j=0;j<N;j++)
				{
					h0(i,j)=std::stod(tokens[j]);
				}
			}
		}
		else
		{
			return h0;
		}
	}
	return h0;
}
