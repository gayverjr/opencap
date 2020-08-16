/*Copyright (c) 2020 James Gayvert

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

#include "Shell.h"
#include <iostream>
#include <vector>
#include <map>
#include <Eigen/Dense>

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

void fill_mat(std::vector<double> &matrix_elements, Eigen::MatrixXd &opdm)
{
	size_t vec_idx = 0;
	for (size_t row_idx=0;row_idx<opdm.rows();row_idx++)
	{
		for (size_t col_idx=0;col_idx<opdm.cols();col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
	}
}

void fill_LT(std::vector<double> matrix_elements, Eigen::MatrixXd &opdm)
{
	size_t vec_idx = 0;
	size_t row_idx = 0;
	while(row_idx<opdm.rows() && vec_idx<matrix_elements.size())
	{
		//elements are added to each column <= row index
		for (size_t col_idx=0;col_idx<=row_idx;col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			opdm(col_idx,row_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
		row_idx++;
	}
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

