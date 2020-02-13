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

std::vector<std::array<size_t,3>> get_carts_ordering(Shell shell)
{
	std::vector<std::array<size_t,3>> ordering;
	size_t a=shell.l; size_t b=0; size_t c=0;
	while (ordering.size()<shell.num_carts())
	{
		ordering.push_back({{a,b,c}});
		if (c<shell.l-a)
		{
			b=b-1;
			c=c+1;
		}
		else
		{
			a=a-1;
			c=0;
			b=shell.l-a;
		}
	}
	return ordering;
}

std::vector<int> get_harmonic_ordering(Shell shell)
{
	std::vector<int> order;
	if (shell.l ==0)
		return {0};
	for(int i=-1*shell.l;i<shell.l;i++)
	{
		order.push_back(i);
	}
	return order;
}


