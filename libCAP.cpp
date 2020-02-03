#include <iostream>
#include "BasisFunction.h"
#include <vector>
#include "overlap.h"
#include "BasisSet.h"
using namespace std;

int main(int argc, char **argv) {

	//cout << "Hello world" << std::endl;
	BasisSet bs("N2.xyz","bas.bas");
	//cout << bs.Nbasis << std::endl;
	std::vector<double> exps1 = {69.12110,15.83500,4.673260};
	std::vector<double> coeffs1 = {0.0356574,0.2394770,0.8184610};
	std::vector<double> exps2={2.304};
	std::vector<double> coeffs2 = {1.0};
	std::vector<int> shell1 = {0,0,0};
	std::vector<int> shell2 = {0,0,2};
	std::vector<double> coords = {0.0,0.0,0.0};
	std::vector<double> coords2 = {1.0,1.0,1.0};
	BasisFunction bf1 = BasisFunction(coords,shell1,exps1,coeffs1);
	BasisFunction bf2 = BasisFunction(coords2,shell2,exps2,coeffs2);
	//cout << overlap_integral(bf1,bf2);

	cout << bs.Nbasis << std::endl;
	cout << "5xx and 5z^2" << std::endl << overlap_integral(bf1,bf2) << std::endl;
	for(int i=0;i<bs.Nbasis;i++)
	{
		for(int j=0;j<bs.Nbasis;j++)
		{
			double res = overlap_integral(bs.basis[i],bs.basis[j]);
			cout << i+1 << " " << j+1 << std::endl;
			cout << res << std::endl;
		}
	}

	return 0;
}


