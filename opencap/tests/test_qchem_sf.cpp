/*Copyright (c) 2020 James Gayvert

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

/*
 * test_qchem_sf.cpp
 *
 */
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


int main()
{
	std::string input_filename = "../tests/qchem/test_qchem_sf.in";
	try
	{
		std::tuple<System,std::map<std::string,std::string>> inp_data = parse_input(input_filename);
		std::map<std::string,std::string> params = std::get<1>(inp_data);
		if(params["jobtype"] == "perturb_cap")
		{
			CAP pc(std::get<0>(inp_data),get_params_for_field(params,"perturb_cap"));
			pc.run();
		}
	}
	catch (exception& e)
	{
		opencap_handle_exception(e);
		return 1;
	}
	return 0;
}


