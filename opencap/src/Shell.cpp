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

#include "Shell.h"

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "Atom.h"
#include "utils.h"

using namespace std;

Shell::Shell(int angmom, bool cart_flag,std::vector<double> exponents,
		std::vector<double> coefficients)
{
	l = angmom;
	pure = cart_flag;
	exps=exponents;
	coeffs=coefficients;
	num_prims = exps.size();
	origin = {{0.0,0.0,0.0}};
	normalize();
}

Shell::Shell(int angmom,std::array<double,3>center)
{
	l = angmom;
	num_prims = 0;
	pure = true;
	origin = center;
}

Shell::Shell()
{
	l = -1;
	num_prims = 0;
	pure = true;
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

size_t Shell::num_bf()
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

void Shell::evaluate_on_grid(double* x, double* y, double* z, int num_points,size_t lx,size_t ly, size_t lz, Eigen::Ref<Eigen::VectorXd> v)
{
    #pragma omp parallel for
	for (size_t i=0;i<num_points;i++)
	{
		double result = 0;
		for(size_t j=0;j<num_prims;j++)
		{
			double r_squared = pow(x[i]-origin[0],2) + pow(y[i]-origin[1],2) + pow(z[i]-origin[2],2);
			result+=coeffs[j]*pow(x[i]-origin[0],lx)*pow(y[i]-origin[1],ly)*pow(z[i]-origin[2],lz)*pow(euler,-1.0*r_squared*exps[j]);
		}
		v(i) = result;
	}

}

void Shell::add_primitive(double exp,double coeff)
{
	exps.push_back(exp);
	coeffs.push_back(coeff);
	num_prims = exps.size();
}

bool Shell::operator==(const Shell& other)
{
	return l == other.l && exps == other.exps && coeffs==other.coeffs && origin==other.origin;
}







