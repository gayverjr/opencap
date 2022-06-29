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

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
//#include <boost/math/special_functions/gamma.hpp>  remove for now, keep in case we ever end up using boost
#include <vector>
#include "Shell.h"
#include "boxcap.h"
#include "../asa/asa239.h"

using namespace std;


const double PI = 3.14159265358979;
const double sqrt_pi = 1.77245385090552;
const double OOPI3 = 1.0 / (PI * PI * PI);
const double gauss_factor[] = {1.0, 0.5, 0.75, 1.875, 6.5625, 29.53125, 162.421875};

//
// box-CAP helper functions
//
double intmod_r2(double a, int k, double cl, double cu);
int xi (double z, int i);
double intgauss(double alpha, int n, int = 0, const char* = "message ");
vector<double> polmulti(vector<double>& g, int n, vector<double>& h, int m,
	       int = 0, const char* = "message");
void poltrans (vector<double>& f, double alpha, double beta, int n,
	       int = 0, const char* = "message ");
void gauss2(double ao, double alpha, double bo, double beta,
	    double& co, double& expo, double& factor, int = 0);
void print_vector(vector<double>& v, const char msg[]);

/*
     gaussian integral of the softbox potential 
     V  = { (x+boxlength)^2  for x < -boxlength
          {        0         for -boxlength < x < boxlength
          { (x-boxlength)^2  for x > boxlength 

This is a Fortran to Cpp translation of our old code from the late 1990ties.
Fortran by Uwe Riss, Robin Santra, and Thomas Sommerfeld
Cpp translation by Thomas Sommerfeld, Nov 2018 T. Sommerfeld

The integrals of the so-called box-CAP are computed analytically.
Do not even try to backwards-engineer this, dig out Robin's paper. 

@Article{santra99,
author = "R. Santra and L. S. Cederbaum and H.-D. Meyer",
journal = "Chem. Phys. Lett.",
volume = 303,
year =   1999,
pages =  413,
title =  "Electronic decay of molecular clusters: Non-stationary states computed by standard quantum chemistry methods"}

gto1 and gto2 define two primitive Cartesian GTOs including their l value,
say, px, or dx^2.
boxlength[3] are the cutoff parameters defined above
normalized [0] is a flag 0=unnormalized GTOs, 1=normalized GTOs
verbose can be increased to create debug information  

*/

double integrate_box_cap(Shell shell1, Shell shell2, 
std::array<size_t,3> l1, std::array<size_t,3> l2, double boxlength[3],double thresh)
{
  double sum = 0;
	for(size_t prim1 =0;prim1<shell1.num_prims;prim1++)
  {
    gto gto1 = {shell1.origin,shell1.exps[prim1],l1};
    for(size_t prim2=0;prim2<shell2.num_prims;prim2++)
    {
      gto gto2 = {shell2.origin,shell2.exps[prim2],l2};
      if (shell1.exps[prim1] > thresh && shell2.exps[prim2] > thresh)
        sum+= shell1.coeffs[prim1] * shell2.coeffs[prim2] * boxcap(gto1,gto2,boxlength,false,0);
    }
  }
  return sum;
}

double boxcap(gto gto1, gto gto2, double boxlength[3],
	      int normalized, int verbose) {

  int debug = verbose;
  
  double cen[3] = {0, 0, 0};
  double alpha = gto1.a;
  double beta  = gto2.a;
  double d0, delta, fakt;
  
  double ovl[3] = {0,0,0};
  double dip[3] = {0,0,0};
  for (int k = 0; k < 3; ++k) {
    int m = gto1.l[k];
    int n = gto2.l[k];
    int nm = n+m;
    double ao = gto1.r[k];
    double bo = gto2.r[k];
    gauss2(ao, alpha,  bo, beta, d0, delta, fakt, verbose);    
    
    vector<double> g (m+1,0);
    vector<double> h (n+1,0);
    g[m] = 1;
    h[n] = 1;
    poltrans(g, ao, d0, m, debug, " g = ");
    poltrans(h, bo, d0, n, debug, " h = ");
    vector<double> f = polmulti(g, m, h, n, debug, " f = ");

    double cl = cen[k] - d0 - boxlength[k];
    double cu = cen[k] - d0 + boxlength[k];
    
    double temp1 = 0;
    double temp2 = 0;
    for (int i = 0; i <= nm; ++i) 
    {
      temp1 += f[i] * intgauss( delta, i, debug, "  I = ");
      temp2 += f[i] * intmod_r2(delta, i, cl, cu);
    }
    ovl[k] = fakt * temp1;
    dip[k] = fakt * temp2;
  }
  if (verbose > 1) {
    cout << "Ovl: " << ovl[0] << ", " <<  ovl[1] << ", " <<  ovl[2] << endl;
    cout << "Dip: " << dip[0] << ", " <<  dip[1] << ", " <<  dip[2] << endl; 
  }

  double cap = dip[0]*ovl[1]*ovl[2] + ovl[0]*dip[1]*ovl[2] + ovl[0]*ovl[1]*dip[2];
  if (normalized == 0)
    return cap;
  else {
    alpha *= 2;
    beta *= 2;
    double ovl[3] = {0, 0, 0};
    for (int k=0; k < 3; ++k) {
      int m = 2 * gto1.l[k]; 
      int n = 2 * gto2.l[k];
      ovl[k] = intgauss(alpha,m) * intgauss(beta,n);
    }
    if (verbose > 1)
      cout << "<0>: " << ovl[0] << ", " <<  ovl[1] << ", " <<  ovl[2] << endl;
    return cap / sqrt(ovl[0]*ovl[1]*ovl[2]);
  }
}


//
// 
// integrate(x**k * f(x) * exp(-a*x**2), x, -oo, oo)
//
// where 
//
//        (x + cl)**2 if x < cl
// f(x) =    0        if cl <= x <= cu 
//        (x - cu)**2 if x > cu
//
// This integral can be expressed in terms of the incomplete gamma-function.    
// see:
//  @Article{santra99,
//  author = "R. Santra and L. S. Cederbaum and H.-D. Meyer",
//  journal = "Chem. Phys. Lett.",
//  volume = 303,
//  year =   1999,
//  pages =  413,
//  title =  "Electronic decay of molecular clusters:
//            Non-stationary states computed by standard quantum chemistry methods"}
//
//
double intmod_r2(double a, int k, double cl, double cu) {
  int sg = -1;
  if ((k % 2) == 0)
    sg = 1;
  double ql = a*cl*cl;
  double qu = a*cu*cu;
  double result = 0;
  int ifault = 0;
  for (int i = 0; i < 3; ++i) {
    int j = k + i;
    double par = 0.5*(j + 1);
    int ex = 2 - i;
    double dl = pow(cl,ex);
    double du = pow(-cu,ex);
    double term = tgamma(par)*(sg*dl + du);
    double norm_factor = tgamma(par);
    term -= sg*dl*gammad(ql,par,&ifault)*norm_factor*xi(-cl,j);
    term -=    du*gammad(qu,par,&ifault)*norm_factor*xi( cu,j);
    //term -= sg*dl*boost::math::tgamma_lower(par,ql)*xi(-cl,j);  func is not normalized w/ boost
    //term -=    du*boost::math::tgamma_lower(par,qu)*xi( cu,j);  func is not normalized w/ boost
    term /= pow(a,par);
    if (i == 1)
      term *= 2;
    result += term;
  }
  return 0.5*result;
}


//
//  computes a sign needed in intmod_r2()
//
int xi (double z, int i) {
  if (z >= 0)
    return 1;
  else {
    if ((i%2) == 0)
      return -1;
    else
      return 1;
  }
}




//
//  returns int dx x^n exp(-alpha*x^2)
//
double intgauss(double alpha, int n, int debug, const char* msg) {
  double integral = 0.0;
  if (n%2 == 0) {
    int k = n/2;
    double a = pow(alpha,k);
    integral = gauss_factor[k] * sqrt(PI / alpha) / a;
  }
  if (debug > 3)
    cout << msg << integral << endl;
  return integral;
}


//
//   {sum} a_i x^i  *  {sum} b_i x^i  ->  {sum} c_i x^i 
//
vector<double> polmulti(vector<double>& g, int n, vector<double>& h, int m,
			int verbose, const char msg[]) {
  vector<double> f (n+m+1,0);
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
      f[i+j] += g[i] * h[j];
    }
  }
  if (verbose > 2) {
    print_vector(f, msg);
  }
  return f;
}


//
// {sum} a_i ( x - a )^i   ->   {sum} b_i ( x - b )^i
//
void poltrans (vector<double> &f, double alpha, double beta, int n,
	       int verbose, const char msg[]) {
  double diff = beta - alpha;
  int k = 1;
  for (int j = 0; j <= n; ++j) {
    double w = 0;
    for (int i = n; i >= j; --i) {
      w = w*diff + f[i];
      f[i] *= (i-j);
    }
    f[j] = w / k;
    k = (j+1)*k;
  }
  if (verbose > 2)
    print_vector(f, msg);
}


void gauss2(double ao, double alpha, double bo, double beta,
	    double& co, double& expo, double& factor, int verbose) {
  expo = alpha + beta;
  co = (ao*alpha + bo*beta) / expo;
  factor = alpha * beta * (bo-ao) * (bo-ao) / expo;
  factor = exp(-factor);
  if (verbose > 2)
    cout << "gs2:" << co << "  " << expo << "  "  << factor << endl;
}


void print_vector(vector<double> &v, const char msg[]) {
  int n = static_cast<int>(v.size());
  cout << msg;
  for (int i = 0; i < n; ++i) {
    cout << v[i];
    if ((i+1)%5 == 0)
      cout << endl;
    else
      cout << ", ";
  }
  if (n%5 != 0)
    cout << endl;
}

