/*Copyright (c) 2021 James Gayvert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "Atom.h"
#include "CAP.h"
#include "date.h"
#include "InputParser.h"
#include "keywords.h"
#include "molcas_interface.h"
#include "opencap_exception.h"
#include "System.h"
#include "TrajectoryAnalysis.h"

int main(int argc, char **argv)
{
	std::cout<< "Welcome to OpenCAP: An open-source program for studying resonances in molecules." << std::endl
			 << "Version 1.1.1" << std::endl
			 << "Developed by James Gayvert and Ksenia Bravaya. Distributed under the MIT license." << std::endl
			 << "For questions and support, please visit: https://github.com/gayverjr/opencap"
			 <<  std::endl;
	if (argc == 2)
	{
		std::string input_filename = argv[1];
		try
		{
			std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input(input_filename);
			std::map<std::string,std::string> params = std::get<1>(inp_data);
			if(params["projected_cap"] == "true")
			{
				std::cout << std::endl << "Starting Projected CAP calculation." << std::endl;
				auto t_start = std::chrono::high_resolution_clock::now();
				CAP pc(std::get<0>(inp_data),get_params_for_field(params,"projected_cap"));
				pc.run();
				std::cout << "Printing out matrices required for Projected CAP calculation." << std::endl;
				std::cout << "Number of states: " << pc.nstates << std::endl;
				std::cout << "Zeroth order Hamiltonian" << std::endl;
				std::cout << std::setprecision(17)<< pc.ZERO_ORDER_H << std::endl;
				std::cout << "CAP matrix" << std::endl;
				std::cout << std::setprecision(17) << std::scientific << pc.CAP_MAT << std::endl;
				auto t_end = std::chrono::high_resolution_clock::now();
				std::cout << std::fixed << "Wall time:" << std::chrono::duration<double>(t_end-t_start).count() << std::endl;
				if(params["trajectory"]=="true")
				{
					CAPHamiltonian caph(pc.ZERO_ORDER_H,pc.CAP_MAT,get_params_for_field(params,"trajectory"),input_filename);
					caph.run_trajectory();
				}
				// print time
			    using namespace date;
			    using namespace std::chrono;
				auto tp = system_clock::now();
				const auto tpm = floor<minutes>(tp);
				const auto dp = floor<days>(tpm);
				const auto ymd = year_month_day{dp};
				auto time = make_time(tpm-dp);
				std::cout << "Job finished: " << ymd << ' ' << time << std::endl;
				std::cout << "Thank you for using OpenCAP!" << std::endl;
			}
		}
		catch (exception& e)
		{
			opencap_handle_exception(e);
		}
	}
	else

		std::cout << "Usage: opencap infile" << std::endl;
	return 0;
}

