/*
 * transforms.cpp
 *
 *  Created on: Feb 12, 2020
 *      Author: JG
 */
#include <armadillo>
#include "BasisSet.h"
#include "Shell.h"
#include "utils.h"
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include "gto_ordering.h"


void uniform_cart_norm(arma::mat &my_mat, BasisSet bs)
{
	unsigned int bf_idx = 0;
	for(auto&shell:bs.basis)
	{
		std::vector<std::array<size_t,3>> order = opencap_carts_ordering(shell);
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

arma::mat get_trans_mat(Shell shell)
{
	std::vector<std::array<size_t,3>> cart_order = opencap_carts_ordering(shell);
	std::vector<int> sph_order = opencap_harmonic_ordering(shell);
	arma::mat trans_mat(shell.num_bf,shell.num_carts());
	for(size_t i=0;i<shell.num_bf;i++)
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

void transform_block(Shell shell1, Shell shell2, arma::subview<double>&cart_sub_mat,arma::subview<double>&sph_sub_mat)
{
	if(!shell1.pure && !shell2.pure)
		sph_sub_mat = cart_sub_mat;
	else if(shell1.pure && !shell2.pure)
		sph_sub_mat = get_trans_mat(shell1)*cart_sub_mat;
	else if(shell1.pure && shell2.pure)
		sph_sub_mat = get_trans_mat(shell1)*cart_sub_mat*get_trans_mat(shell2).t();
	else
		sph_sub_mat = cart_sub_mat*get_trans_mat(shell2).t();
}

void cart2spherical(arma::mat &cart_ints, arma::mat &spherical_ints, BasisSet bs)
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
			auto cart_sub_mat = cart_ints.submat(cart_row_idx,cart_col_idx,
								cart_row_idx+shell1.num_carts()-1,cart_col_idx+shell2.num_carts()-1);
			auto sph_sub_mat = spherical_ints.submat(sph_row_idx,sph_col_idx,
					sph_row_idx+shell1.num_bf-1,sph_col_idx+shell2.num_bf-1);
			transform_block(shell1,shell2,cart_sub_mat,sph_sub_mat);
			cart_col_idx+=shell2.num_carts();
			sph_col_idx+=shell2.num_bf;
		}
		cart_row_idx+=shell1.num_carts();
		sph_row_idx+=shell1.num_bf;
	}
}

