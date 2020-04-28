#include <iostream>
#include <vector>
#include "System.h"
#include <chrono>
#include <map>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>
#include "InputParser.h"
#include "Atom.h"

int main(int argc, char **argv)
{
	std::cout << "Welcome to OpenCAP!" << std::endl;
	auto t_start = std::chrono::high_resolution_clock::now();
	std::tuple<std::vector<Atom>,std::map<std::string,std::string>> inp_data = parse_input(argv[1]);
	std::map<std::string,std::string> params = std::get<1>(inp_data);
	System my_sys(std::get<0>(inp_data),std::get<1>(inp_data));
	my_sys.compute_cap_matrix();
	std::cout << "Printing out matrices required for Projected CAP calculation:" << std::endl;
	std::cout << "Number of states: " << my_sys.nstates << std::endl;
	std::cout << "Zeroth order Hamiltonian" << std::endl;
	my_sys.ZERO_ORDER_H.raw_print();
	std::cout << "CAP matrix" << std::endl;
	my_sys.CORRELATED_CAP_MAT.raw_print();
	auto t_end = std::chrono::high_resolution_clock::now();
	std::cout << "Wall time:" << std::chrono::duration<double, std::milli>(t_end-t_start).count() << std::endl;
	return 0;
}

