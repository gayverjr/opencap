#include <iostream>
#include <vector>
#include "BasisFunction.h"
#include <math.h>
# define M_PIl          3.141592653589793238462643383279502884
using namespace std;

BasisFunction::BasisFunction(std::vector<double> coords, std::vector<int> shell_vec,
		std::vector<double> basis_exps, std::vector<double> basis_coeffs)
{
	center = coords;
	shell = shell_vec;
	exps = basis_exps;
	coeffs= basis_coeffs;
	num_prims = exps.size();
	normalize();
}

int BasisFunction::fact2(unsigned int n)
{
    int res = 1;
    for (int i=n; i>=0; i=i-2)
    {
        if (i==0 || i==1)
            return res;
        else
            res *= i;
    }
    return res;
}

void BasisFunction::normalize()
{	//G_ijk = x^l * y^m * z^n * exp(-ar^2)
	int L = shell[0] + shell[1] + shell[2]; //total angular momentum quantum number
	//normalize primitives first
	for(unsigned int i =0;i<exps.size();i++)
	{
		double res = sqrt(pow(4*exps[i],L)*pow(2*exps[i]/M_PIl,1.5)
				/fact2(2*shell[0]-1)/fact2(2*shell[1]-1)/fact2(2*shell[2]-1));
		norms.push_back(res);
	}
	//pi^(3/2)*(2l-1)!!*(2m-1)!!*(2n-1)!!/2^L
	double prefactor = pow(M_PIl,1.5)*fact2(2*shell[0]-1)*fact2(2*shell[1]-1)*fact2(2*shell[2]-1)/pow(2,L);
	double N = 0;
	for(int i=0;i<num_prims;i++)
	{
		for(int j=0;j<num_prims;j++)
		{
			N+=norms[i]*coeffs[i]*norms[j]*coeffs[j]/pow(exps[i]+exps[j],L+1.5);
		}
	}
	N*=prefactor;
	N=pow(N,-0.5);
	for(int i=0; i<num_prims;i++)
	{
		coeffs[i]*=N;
	}

}
