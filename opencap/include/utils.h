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

/*! \file utils.h
     \brief Helper functions used by various pieces of the software.
 */
#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <map>


//useful constants
#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884
#endif
#define euler 2.718281828459045235360287471352662497
/*! Computes double factorial (https://mathworld.wolfram.com/DoubleFactorial.html)
 *  \param n : integer to be evaluated
*/
double fact2(int n);
/*! Computes factorial
 * \param n : integer to be evaluated
*/
double fact(int n);
/*! Returns -1 if odd, 1 if even
 *  \param i : integer to be evaluated
*/
int parity(int i);
/*! Computes binomial coefficient for P choose Q
 *  \param p : integer to be evaluated
 *  \param q : integer to be evaluated
*/
double binom(int p, int q);
/*! Helper function for splitting strings
    \param s: string to be split,
    \param delimiter: delimiter for splitting
    \return A vector containing the split string
*/
std::vector<std::string> split(const std::string& s, char delimiter);
/*! Populates symmetric matrix given list of elements read in row major order in lower triangular form
 * \param matrix_elements: list of lower triangular matrix elements in row-major order
 * \param opdm: matrix to be populated
*/

/*
void fill_LT(std::vector<double> matrix_elements, Eigen::MatrixXd &opdm);
/*! Populates matrix given full list of elements in row major order
 * \param matrix_elements: list of all matrix elements in row-major order
 * \param opdm: matrix to be populated
*/
/*
void fill_mat(std::vector<double> &matrix_elements, Eigen::MatrixXd &opdm);
/*! Converts fortran-style d floats to e c-style floats.
 */


template <typename T>
void fill_mat(std::vector<T> &matrix_elements, Eigen::MatrixXd &opdm)
{
	size_t vec_idx = 0;
	for (size_t row_idx=0;row_idx<opdm.rows();row_idx++)
	{
		for (size_t col_idx=0;col_idx<opdm.cols();col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
	}
}

template <typename T>
void fill_LT(std::vector<T> matrix_elements, Eigen::MatrixXd &opdm)
{
	size_t vec_idx = 0;
	size_t row_idx = 0;
	while(row_idx<opdm.rows() && vec_idx<matrix_elements.size())
	{
		//elements are added to each column <= row index
		for (size_t col_idx=0;col_idx<=row_idx;col_idx++)
		{
			opdm(row_idx,col_idx) = matrix_elements[vec_idx];
			opdm(col_idx,row_idx) = matrix_elements[vec_idx];
			vec_idx++;
		}
		row_idx++;
	}
}

void fortran_dfloats_to_efloats(std::string &str);
/*! Converts angular momentum letter(s,p,d etc.) to quantum number l
 */
int shell2angmom(std::string angmom);
/*! Case insensitive comparison of two strings
 */
bool compare_strings(std::string s1, std::string s2);
/*! Checks if the string is an integer
 */
bool is_number(const std::string &s);
/*! Checks if the string is a letter
 */
bool is_letter(const std::string &s);
/*! Reads formatted matrix from .txt file
 */
Eigen::MatrixXd read_matrix(size_t N, std::string filename);

#endif /* UTILS_H_ */
