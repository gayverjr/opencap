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

#include "overlap.h"

#include <math.h>
#include <cstdlib>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "gto_ordering.h"
#include "utils.h"

//overlap between pair of gaussians
double gauss_integral(double exp_a,std::array<size_t, 3> cart_a,std::array<double,3> coord_a,
		double exp_b,std::array<size_t, 3> cart_b, std::array<double,3> coord_b)
{
	std::vector<double> integrals;
	for(int i=0;i<3;i++)
	{
		double dist = coord_a[i]-coord_b[i];
		double res = mcmurchie_davidson(exp_a,exp_b,0,cart_a[i],cart_b[i],dist);
		integrals.push_back(res);
	}
	double int_tot = 1.0;
	for (unsigned int i=0;i<integrals.size();i++)
	{
		int_tot*=integrals[i];
	}
	return int_tot*pow(M_PIl/(exp_a+exp_b),1.5);
}

//recursive formula for cartesian overlap distributions
double mcmurchie_davidson(double exp_a, double exp_b, int t, int angmom_a, int angmom_b, double dist)
{
	double p = exp_a + exp_b;
	double q = exp_a*exp_b/p;
	if (t<0 || t>(angmom_a+angmom_b))
	{
		return 0.0;
	}
	else if(angmom_a==angmom_b && angmom_b==0)
	{
		return exp(-q*dist*dist);
	}
	else if(angmom_b==0)
	{
		//decrement angmom_a
		return(
				(1/(2*p))*mcmurchie_davidson(exp_a,exp_b,t-1,angmom_a-1,angmom_b,dist)   -
				(q*dist/exp_a)*mcmurchie_davidson(exp_a,exp_b,t,angmom_a-1,angmom_b,dist) +
				(t+1)*mcmurchie_davidson(exp_a,exp_b,t+1,angmom_a-1,angmom_b,dist));
	}
	else
	{
		//decrement angmom_b
		return(
				(1/(2*p))*mcmurchie_davidson(exp_a,exp_b,t-1,angmom_a,angmom_b-1,dist)+
				(q*dist/exp_b)*mcmurchie_davidson(exp_a,exp_b,t,angmom_a,angmom_b-1,dist)+
				(t+1)*mcmurchie_davidson(exp_a,exp_b,t+1,angmom_a,angmom_b-1,dist));
	}
}

//overlap between two normalized, contracted gaussians
double overlap_integral(Shell a, std::array<size_t,3> cart_a, Shell b,
		std::array<size_t,3> cart_b)
{
	double res = 0.0;
	for(size_t i=0;i<a.num_prims;i++)
	{
		for(size_t j=0;j<b.num_prims;j++)
		{
			res+=a.coeffs[i]*b.coeffs[j]*
			gauss_integral(a.exps[i],cart_a,a.origin,b.exps[j],cart_b,b.origin);
		}
	}
	return res;
}

Eigen::MatrixXd shell_overlap(Shell shell_a, Shell shell_b)
{
	Eigen::MatrixXd block(shell_a.num_carts(),shell_b.num_carts());
	std::vector<std::array<size_t,3>> order_a = opencap_carts_ordering(shell_a.l);
	std::vector<std::array<size_t,3>> order_b = opencap_carts_ordering(shell_b.l);
	for(size_t i=0;i<shell_a.num_carts();i++)
	{
		std::array<size_t,3> a_cart = order_a[i];
		for(size_t j=0;j<shell_b.num_carts();j++)
		{
			std::array<size_t,3> b_cart = order_b[j];
			block(i,j) = overlap_integral(shell_a,a_cart, shell_b, b_cart);
		}
	}
	return block;
}

void compute_analytical_overlap(BasisSet bs, Eigen::MatrixXd &Smat)
{
	unsigned int row_idx = 0;
	for(auto&shell1:bs.basis)
	{
		unsigned int col_idx = 0;
		for(auto&shell2: bs.basis)
		{
			Smat.block(row_idx,col_idx,shell1.num_carts(),shell2.num_carts()) = shell_overlap(shell1,shell2);
            col_idx += shell2.num_carts();
		}
		row_idx += shell1.num_carts();
	}
}



