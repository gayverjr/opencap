#include "Shell.h"
#include <iostream>
#include <vector>
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
	else
		return {0,0,0};
}




