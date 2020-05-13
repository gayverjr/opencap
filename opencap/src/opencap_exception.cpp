#include "opencap_exception.h"
#include <iostream>
#include <string>
#include <cstdlib>

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
