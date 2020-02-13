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
#include <cmath>


void uniform_cart_norm(arma::mat &my_mat, BasisSet bs)
{
	unsigned int bf_idx = 0;
	for(auto&shell:bs.basis)
	{
		std::vector<std::array<size_t,3>> order = get_carts_ordering(shell);
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

//for m<0, we want imaginary coefficients
//expression is 1j^(n)
//so if m <0 , we return zero if n%2=0, negative if n+1%4=0, positive if n+1%4!=0,
//if m>=0, we return zero if n%2!=0, positive if n%4=0, negative otherwise
double real_imag_coeff_parity(int m, int n)
{
	if(n==0)
		return 1;
	if(m<0)
	{
		if(n%2==0)
			return 0;
		else if((n+1)%4==0)
			return -1;
		else
			return 1;
	}
	else
	{
		if(n%2!=0)
			return 0;
		else if(n%4==0)
			return 1;
		else
			return -1;
	}
}

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
      auto term3=0;
      for (int k=0;k<=j;k++)
      {
    	  term3+=binom(j,k)*binom(abs_m,lx-2*k)
    			  *real_imag_coeff_parity(m,(abs_m - lx +2*k)/2);
      }
      return term1*term2*term3;
}

void transform_block(arma::subview<double>&row,unsigned int L,
		unsigned int m, std::vector<std::array<size_t,3>> cart_order)
{
	for(unsigned int i=0;i<cart_order.size();i++)
	{
		auto cart_angmom = cart_order[i];
		row(0,i) = get_coeff(L,m,cart_angmom[0],cart_angmom[1],cart_angmom[2]);
	}
}

void cart2spherical(arma::mat &cart_ints, BasisSet bs)
{
	arma::mat trans_mat(bs.Nbasis,bs.num_carts());
	trans_mat.zeros();
	unsigned int bf_idx = 0;
	for(auto shell:bs.basis)
	{
		std::vector<std::array<size_t,3>> cart_order = get_carts_ordering(shell);
		std::vector<int> harmonic_order = get_harmonic_ordering(shell);
		for(size_t i=0;i<shell.num_bf;i++)
		{
			auto row = trans_mat.submat(bf_idx+i,0,bf_idx+i,shell.num_carts()-1);
			transform_block(row,shell.l,harmonic_order[i],cart_order);
		}
	}
	arma::mat new_mat = trans_mat * cart_ints;
	new_mat.print();
	std::cout << std::endl;

}


