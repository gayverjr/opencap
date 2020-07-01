#include "Shell.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
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

bool same_atom(std::array<double,3> shell_origin,std::array<double,3> atom_coords)
{
	for(size_t i=0;i<3;i++)
	{
		if(shell_origin[i]!=atom_coords[i])
			return false;
	}
	return true;
}


