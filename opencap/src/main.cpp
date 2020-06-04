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
#include "ProjectedCAP.h"
#include "keywords.h"

int main(int argc, char **argv)
{
	std::cout << "Welcome to OpenCAP!" << std::endl;
	if (argc == 2)
	{
		std::string input_filename = argv[1];
		try
		{
			std::tuple<std::vector<Atom>,std::map<std::string,std::string>> inp_data = parse_input(input_filename);
			std::map<std::string,std::string> params = std::get<1>(inp_data);
			System my_sys(std::get<0>(inp_data),get_params_for_field(params,"system"));
			if(params["jobtype"] == "projected_cap")
			{
				auto t_start = std::chrono::high_resolution_clock::now();
				Projected_CAP pc(my_sys,get_params_for_field(params,"projected_cap"));
				pc.run();
				std::cout << "Printing out matrices required for Projected CAP calculation." << std::endl;
				std::cout << "Number of states: " << pc.nstates << std::endl;
				std::cout << "Zeroth order Hamiltonian" << std::endl;
				pc.ZERO_ORDER_H.raw_print();
				std::cout << "CAP matrix" << std::endl;
				pc.CORRELATED_CAP_MAT.raw_print();
				auto t_end = std::chrono::high_resolution_clock::now();
				std::cout << "Wall time:" << std::chrono::duration<double>(t_end-t_start).count() << std::endl;
			}
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

