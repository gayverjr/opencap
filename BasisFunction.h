#include <vector>
#pragma once
using namespace std;

class BasisFunction {
public:
    std::vector<double>center;
    std::vector<int> shell;
    std::vector<double>exps;
    std::vector<double> coeffs;
    std::vector<double> norms;
    int num_prims;

  public:
	BasisFunction(std::vector<double> coords, std::vector<int> shell_vec,
		std::vector<double> basis_exps, std::vector<double> basis_coeffs);
  	int fact2(unsigned int n);
  private:
  	void normalize();

};
