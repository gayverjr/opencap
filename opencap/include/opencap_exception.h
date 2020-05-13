 /*! \file opencap_exception.h
     Exception handling adapted from: https://github.com/GPMueller/mwe-cpp-exception
 */
#pragma once
#include <stdexcept>
#include <string>

/*! \brief Class for handling nested exceptions.
 *
 */
class Exception : public std::runtime_error
{
		public:
        Exception(const std::string & message, const char * file, unsigned int line):
            std::runtime_error(message)
        {
            _message = std::string(file) + ":" + std::to_string(line) + " : " + message;
        }

        ~Exception() throw() {}

        const char * what() const throw()
        {
            return _message.c_str();
        }

    private:
        std::string _message;
 };

/** Rethrow (creates a std::nested_exception) an exception, using the Exception class
 * which contains file and line info. The original exception is preserved...
 */
void rethrow(const std::string & message, const char * file, unsigned int line);

/** General exception handler
 */
void Handle_Exception(const std::exception & ex, const std::string & function);

// Shorthand for throwing an Exception with file and line info using macros
#define opencap_throw(message) throw Exception(message, __FILE__, __LINE__);

// Shorthand for rethrowing and Exception with file and line info using macros
#define opencap_rethrow(message) rethrow(message, __FILE__, __LINE__);

// Shorthand for handling an exception, including a backtrace
#define opencap_handle_exception(ex) Handle_Exception(ex, __func__);
