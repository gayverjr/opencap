#include <iostream>
#include <vector>
#include "BasisFunction.h"
#include <math.h>
#include "overlap.h"
#include <cstdlib>
# define M_PIl          3.141592653589793238462643383279502884

//overlap between pair of gaussians
double gauss_integral(double exp_a,std::vector<int> shell_a,std::vector<double> coord_a,
		double exp_b, std::vector<int>shell_b, std::vector<double> coord_b)
{
	std::vector<double> integrals;
	for(int i=0;i<3;i++)
	{
		double dist = coord_a[i]-coord_b[i];
		double res = mcmurchie_davidson(exp_a,exp_b,0,shell_a[i],shell_b[i],dist);
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
	//cout << "P:" << p << std::endl << "Q:" << q << std::endl;
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
				(t+1)*mcmurchie_davidson(exp_a,exp_b,t+1,angmom_a-1,angmom_b,dist)
			  );
	}
	else
	{
		//decrement angmom_b
		return(
				(1/(2*p))*mcmurchie_davidson(exp_a,exp_b,t-1,angmom_a,angmom_b-1,dist)+
				(q*dist/exp_b)*mcmurchie_davidson(exp_a,exp_b,t,angmom_a,angmom_b-1,dist)+
				(t+1)*mcmurchie_davidson(exp_a,exp_b,t+1,angmom_a,angmom_b-1,dist)
	          );
	}

}

//overlap between normalized, contracted gaussians
double overlap_integral(BasisFunction a, BasisFunction b)
{
	double res = 0.0;
	for(int i=0;i<a.num_prims;i++)
	{
		for(int j=0;j<b.num_prims;j++)
		{
			res+=a.norms[i]*a.coeffs[i]*b.norms[j]*b.coeffs[j]*
			gauss_integral(a.exps[i],a.shell,a.center,b.exps[j],b.shell,b.center);
		}
	}
	return res;

}

