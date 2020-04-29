#include <iostream>
#include <vector>
#include <math.h>
#include "Atom.h"
#include "Shell.h"
#include "utils.h"
# define M_PIl          3.141592653589793238462643383279502884
# define euler 			2.718281828459045235360287471352662497
using namespace std;

Shell::Shell(int angmom, bool cart_flag,std::vector<double> exponents,
		std::vector<double> coefficients)
{
	l = angmom;
	pure = cart_flag;
	exps=exponents;
	coeffs=coefficients;
	num_prims = exps.size();
	num_bf = get_size();
	origin = {{0.0,0.0,0.0}};
	normalize();
}

Shell::Shell()
{
	l = 0;
	num_prims = 0;
	num_bf = 0;
	pure = false;
	origin = {{0.0,0.0,0.0}};
}

void Shell::update_coords(std::array<double,3> center)
{
	origin = center;
}

size_t Shell::num_carts()
{
  return (l + 1) * (l + 2) / 2;
}

size_t Shell::get_size()
{
	return pure ? (2 * l + 1) : num_carts();
}

void Shell::normalize()
{
	//normalize primitives first
	for(size_t i =0;i<num_prims;i++)
	{
		double res = sqrt(pow(4*exps[i],l)*pow(2*exps[i]/M_PIl,1.5)
				/fact2(2*l-1));
		coeffs[i]*=res;
	}
	//pi^(3/2)*(2l-1)!!
	double prefactor = pow(M_PIl,1.5)*fact2(2*l-1)/pow(2,l);
	double N = 0;
	for(size_t i=0;i<num_prims;i++)
	{
		for(size_t j=0;j<num_prims;j++)
		{
			N+=coeffs[i]*coeffs[j]/pow(exps[i]+exps[j],l+1.5);
		}
	}
	N*=prefactor;
	N=pow(N,-0.5);
	for(size_t i=0; i<num_prims;i++)
	{
		coeffs[i]*=N;
	}
}

double Shell::evaluate(double x, double y, double z, size_t lx, size_t ly, size_t lz)
{
	double result = 0;
	for(size_t i=0;i<num_prims;i++)
	{
		double r_squared = pow(x-origin[0],2) + pow(y-origin[1],2) + pow(z-origin[2],2);
		result+=coeffs[i]*pow(x-origin[0],lx)*pow(y-origin[1],ly)*pow(z-origin[2],lz)*pow(euler,-1.0*r_squared*exps[i]);
	}
	return result;

}

bool Shell::operator==(const Shell& other)
{
	return l == other.l && exps == other.exps && coeffs==other.coeffs && origin==other.origin;
}





