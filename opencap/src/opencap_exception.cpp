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

// Exception handling adapted from: https://github.com/GPMueller/mwe-cpp-exception.
// Copyright (c) 2016 Gideon Müller

#include "opencap_exception.h"

#include <cstdlib>
#include <iostream>
#include <string>

void rethrow(const std::string & message, const char * file, unsigned int line)
{
		try
		{
			std::rethrow_exception(std::current_exception());
		}
		catch (...)
		{
			std::throw_with_nested(Exception(message, file, line));
		}
}

// Backtrace an exception by recursively unwrapping the nested exceptions
void Backtrace(const std::exception & ex)
{
	try
	{
		std::cerr << ex.what() << std::endl;
		rethrow_if_nested(ex);
	}
	catch( const std::exception & nested_ex )
	{
		Backtrace(nested_ex);
	}
}


// Rethrow (creates a std::nested_exception) an exception, using the Exception class
// which contains file and line info. The original exception is preserved...
 // General Exception handler
void Handle_Exception(const std::exception & ex, const std::string & function)
{
	try
	{
		if (function != "")
			std::cerr << "Exception caught in function \'" << function << "\'" << std::endl;
		std::cerr << "Backtrace:" << std::endl;
		Backtrace(ex);
	}
	catch( ... )
	{
		std::cerr << "Something went super-wrong! TERMINATING!" << std::endl;
		std::exit( EXIT_FAILURE );
	}
}
