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
#include "opencap_exception.h"

int main(int argc, char **argv)
{
	std::cout << "Welcome to OpenCAP!" << std::endl;
	auto t_start = std::chrono::high_resolution_clock::now();
	if (argc == 2)
	{
		std::string input_filename = argv[1];
		std::tuple<std::vector<Atom>,std::map<std::string,std::string>> inp_data = parse_input(input_filename);
		std::map<std::string,std::string> params = std::get<1>(inp_data);
		try
		{
			System my_sys(std::get<0>(inp_data),std::get<1>(inp_data));
			my_sys.compute_cap_matrix();
			std::cout << "Printing out matrices required for Projected CAP calculation." << std::endl;
			std::cout << "Number of states: " << my_sys.nstates << std::endl;
			std::cout << "Zeroth order Hamiltonian" << std::endl;
			my_sys.ZERO_ORDER_H.raw_print();
			std::cout << "CAP matrix" << std::endl;
			my_sys.CORRELATED_CAP_MAT.raw_print();
			auto t_end = std::chrono::high_resolution_clock::now();
			std::cout << "Wall time:" << std::chrono::duration<double>(t_end-t_start).count() << std::endl;
		}
		catch (exception& e)
		{
			opencap_handle_exception(e);
		}
	}
	else
		std::cout << "Error: please specify one input file!" << std::endl;
	return 0;
}

