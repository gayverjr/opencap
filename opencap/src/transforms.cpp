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

/*
 * transforms.cpp
 *
 */

#include <math.h>
#include <algorithm>
#include <array>
#include <Eigen/Dense>
#include <vector>

#include "BasisSet.h"
#include "gto_ordering.h"
#include "Shell.h"
#include "utils.h"

void uniform_cart_norm(Eigen::MatrixXd &my_mat, BasisSet &bs)
{
	unsigned int bf_idx = 0;
	for(auto&shell:bs.basis)
	{
		std::vector<std::array<size_t,3>> order = opencap_carts_ordering(shell.l);
		for(unsigned int i=0;i<shell.num_carts();i++)
		{
			std::array<size_t,3> cart = order[i];
			double scale = sqrt(fact2(2*shell.l-1)/fact2(2*cart[0]-1)
					/fact2(2*cart[1]-1)/fact2(2*cart[2]-1));
			my_mat.col(bf_idx+i) = my_mat.col(bf_idx+i)*scale;
			my_mat.row(bf_idx+i) = my_mat.row(bf_idx+i)*scale;
		}
		bf_idx += shell.num_carts();
	}
}

//DOI: 10.1002/qua.560540202
//Last part of eqn 15, for M<0 we want imaginary part, for M>0 we want real part
double term4(int M, int exp_num)
{
	if(M<0)
	{
		//integer powers will be real
		if(exp_num%2==0)
			return 0;
		else
			return (exp_num-1)%4==0 ? 1:-1;
	}
	else
	{
		//half-integer powers will be imaginary
		if(exp_num%2!=0)
			return 0;
		else
			return exp_num%4==0 ? 1:-1;
	}
}

//DOI: 10.1002/qua.560540202
//Equation 15, currently assumes fully normalized cartesians. Will add non-fully normalized carts soon enough...
double get_coeff(int L, int m, int lx, int ly, int lz)
{
      auto abs_m = std::abs(m);
      if ((lx + ly - abs_m)%2)
        return 0.0;
      auto j = (lx + ly - abs_m)/2;
      if (j < 0)
        return 0.0;
      auto term1 = sqrt( fact(2*lx)*fact(2*ly)*fact(2*lz)*fact(L)*fact(L-abs_m)
    		  	  	  	  /fact(2*L)/fact(lx)/fact(ly)/fact(lz)/fact(L+abs_m));
      term1/=fact(L);
      term1/=pow(2,L);
      auto term2 = 0;
      for(int i=0;i<=(L-abs_m)/2;i++)
      {
    	  term2+=binom(L,i)*binom(i,j)
    			  *parity(i) * fact(2*(L-i))
				  /fact(L-abs_m-2*i);
      }
      double term3=0;
      for (int k=0;k<=j;k++)
    	  term3+=binom(j,k)*binom(abs_m,lx-2*k)*term4(m,abs_m - lx +2*k);
      // for m!=0, real solid harmonics are linear combinations of complex ones
      // R+(l,m) = ( Y(l,m) + Y(l,-m) )/ sqrt(2) ;  R-(l,m) = ( Y(l,m) - Y(l,-m) )/ sqrt(-2)
      double result = (m == 0) ? term1*term2*term3 : M_SQRT2*term1*term2*term3;
      return result;
}

Eigen::MatrixXd get_trans_mat(Shell &shell)
{
	std::vector<std::array<size_t,3>> cart_order = opencap_carts_ordering(shell.l);
	std::vector<int> sph_order = opencap_harmonic_ordering(shell.l);
	Eigen::MatrixXd trans_mat(shell.num_bf(),shell.num_carts());
	for(size_t i=0;i<shell.num_bf();i++)
	{
		int M = sph_order[i];
		for(size_t j=0;j<shell.num_carts();j++)
		{
			std::array<size_t,3> cart = cart_order[j];
			trans_mat(i,j) = get_coeff(shell.l,M,cart[0],cart[1],cart[2]);
		}
	}
	return trans_mat;
}

Eigen::MatrixXd transform_block(Shell &shell1, Shell &shell2, Eigen::MatrixXd cart_block)
{

	if(!shell1.pure && !shell2.pure)
		return cart_block;
	else if(shell1.pure && !shell2.pure)
		return get_trans_mat(shell1)*cart_block;
	else if(shell1.pure && shell2.pure)
		return get_trans_mat(shell1)*cart_block*get_trans_mat(shell2).transpose();
	else
		return cart_block*get_trans_mat(shell2).transpose();
}

void cart2spherical(Eigen::MatrixXd &cart_ints, Eigen::MatrixXd &spherical_ints, BasisSet &bs)
{
	//indices for first basis function for cart and spherical matrices
	unsigned int cart_row_idx = 0;
	unsigned int sph_row_idx = 0;
	for(auto shell1:bs.basis)
	{
		//indices for 2nd basis function for cart and spherical matrices
		unsigned int cart_col_idx = 0;
		unsigned int sph_col_idx = 0;
		for(auto shell2:bs.basis)
		{
			Eigen::MatrixXd cart_block = cart_ints.block(cart_row_idx,cart_col_idx, shell1.num_carts(),shell2.num_carts());
			spherical_ints.block(sph_row_idx,sph_col_idx,shell1.num_bf(),shell2.num_bf())
							= transform_block(shell1,shell2,cart_block);
			cart_col_idx+=shell2.num_carts();
			sph_col_idx+=shell2.num_bf();
		}
		cart_row_idx+=shell1.num_carts();
		sph_row_idx+=shell1.num_bf();
	}
}

