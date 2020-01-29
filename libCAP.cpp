#include <iostream>
#include "BasisFunction.h"
#include <vector>
#include "overlap.h"
#include "BasisSet.h"
using namespace std;

int main(int argc, char **argv) {

	cout << "Hello world" << std::endl;
	std::vector<double> coords = {0.0, 0.0, 0.0 };
	std::vector<int> shell = {1,0,0};
	std::vector<double> exps2= {2.9642,0.803,1.459,1,1};
	std::vector<double> exps = {20.9642,4.803,1.459};
	std::vector<double> coeffs ={0.114660,0.919999,-0.00303068};
	std::vector<double> coeffs2 = {0.114660,0.919999,-0.00303068,0,0};
	BasisFunction bf(coords,shell,exps,coeffs);
	BasisFunction bf2(coords,shell,exps2,coeffs2);
	BasisFunction bf3(coords,shell,exps2,coeffs2);
	cout << overlap_integral(bf3,bf2) << std::endl;
	BasisSet bs("N2.xyz","bas.bas");
	return 0;
}


