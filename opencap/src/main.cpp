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
#include "CAP.h"
#include "keywords.h"
#include "molcas_interface.h"

int main(int argc, char **argv)
{
	std::cout << "Welcome to OpenCAP!" << std::endl;
	if (argc == 2)
	{
		std::string input_filename = argv[1];
		try
		{
			std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input(input_filename);
			std::map<std::string,std::string> params = std::get<1>(inp_data);
			if(params["jobtype"] == "perturb_cap")
			{
				auto t_start = std::chrono::high_resolution_clock::now();
				CAP pc(std::get<0>(inp_data),get_params_for_field(params,"perturb_cap"));
				pc.run();
				std::cout << "Printing out matrices required for Perturbative CAP calculation." << std::endl;
				std::cout << "Number of states: " << pc.nstates << std::endl;
				std::cout << "Zeroth order Hamiltonian" << std::endl;
				std::cout << pc.ZERO_ORDER_H << std::endl;
				std::cout << "CAP matrix" << std::endl;
				std::cout << std::setprecision(8) << std::scientific << pc.CAP_MAT << std::endl;
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

