#include "Shell.h"
#include <iostream>
#include <vector>
#include <armadillo>
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

std::vector<std::array<size_t,3>> libcap_carts_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {{0,0,0}};
	//p
	if(shell.l==1)
		return {{1,0,0},{0,1,0},{0,0,1}};
	//d
	else if (shell.l==2)
		return {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
	//f
	else if (shell.l==3)
		return {{3,0,0},{2,1,0},{2,0,1},{1,0,2},{1,1,1},{1,2,0},{0,3,0},
				{0,2,1},{0,1,2},{0,0,3}};
	//g
	else if (shell.l==4)
		return { {4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},
		         {1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};
	else
	{
		return {{0,0,0}};
	}
}

std::vector<int> libcap_harmonic_ordering(Shell shell)
{
	//s
	if(shell.l == 0)
		return {0};
	//p
	if(shell.l==1)
		return {-1,0,1};
	//d
	else if (shell.l==2)
		return {-2,-1,0,1,2};
	//f
	else if (shell.l==3)
		return {-3,-2,-1,0,1,2,3};
	//g
	else if (shell.l==4)
		return {-4,-3,-2,-1,0,1,2,3,4};
	//h
	else if(shell.l==5)
		return {-5,-4,-3,-2,-1,0,1,2,3,4,5};
	else
		return {0,0,0};
}

void fill_mat(std::vector<double> &matrix_elements, arma::mat &opdm)
{
	size_t vec_idx = 0;
	for (size_t row_idx=0;row_idx<opdm.n_rows;row_idx++)
	{
		for (size_t col_idx=0;col_idx<opdm.n_cols;col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
	}
}

void fill_LT(std::vector<double> matrix_elements, arma::mat &opdm)
{
	size_t vec_idx = 0;
	size_t row_idx = 0;
	while(row_idx<opdm.n_rows && vec_idx<matrix_elements.size())
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



