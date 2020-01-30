#include <vector>
#pragma once
using namespace std;

class Shell {
public:
	string shell;
	std::vector<double> exps;
	std::vector<double> coeffs;
	int num_prims;

public:
	Shell(string shell_label,std::vector<double> exponents, std::vector<double> coefficients);
};
