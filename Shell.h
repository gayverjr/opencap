#include <vector>
#include <array>
#pragma once
using namespace std;

class Shell {
public:
	int l;
	bool pure;
	std::vector<double> exps;
	std::vector<double> coeffs;
	size_t num_bf;
	size_t num_prims;
	std::array<double,3> origin;
	size_t num_carts();
	Shell(int angmom, bool cart_flag,std::vector<double> exponents,
			std::vector<double> coefficients);
	Shell(Shell s, std::array<double,3>center);
	void update_coords(std::array<double,3> center);
	double evaluate(double x, double y, double z,size_t lx,size_t ly, size_t lz);
	std::vector<double> alpha_min();
	double alpha_max();

private:
	void normalize();
	size_t get_size();


};
