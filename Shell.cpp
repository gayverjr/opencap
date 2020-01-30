#include <iostream>
#include <vector>
#include <math.h>
#include "Atom.h"
#include "Shell.h"
using namespace std;

Shell::Shell(string shell_label,std::vector<double> exponents, std::vector<double> coefficients)
{
	shell=shell_label;
	exps=exponents;
	coeffs=coefficients;
	num_prims = coeffs.size();
}
