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

/*! \file opencap_exception.h
     \brief %Exception handling adapted from: https://github.com/GPMueller/mwe-cpp-exception.
     Copyright (c) 2016 Gideon Müller
 */
#pragma once

#include <stdexcept>
#include <string>

/*! \brief Class for handling nested exceptions.
 *  Adapted from: https://github.com/GPMueller/mwe-cpp-exception
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
