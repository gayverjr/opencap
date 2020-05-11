 /*! \file utils.h
     \brief Helper functions used by various pieces of the software.
 */
#include "Shell.h"
#include <armadillo>
/*! Computes double factorial (https://mathworld.wolfram.com/DoubleFactorial.html)
*/
double fact2(int n);
/*! Computes factorial
*/
double fact(int n);
/*! Returns -1 if odd, 1 if even
*/
int parity(int i);
/*! Computes binomial coefficient for P choose Q
*/
double binom(int p, int q);
/*! Helper function for splitting strings
    \param std::string s: string to be split,
    \param char: delimiter for splitting
    \return A vector containing the split string
*/
std::vector<std::string> split(const std::string& s, char delimiter);
/*! Populates symmetric matrix given list of elements read in row major order in lower triangular form
*/
void fill_LT(std::vector<double> matrix_elements, arma::mat &opdm);
/*! Populates matrix given full list of elements in row major order
*/
void fill_mat(std::vector<double> &matrix_elements, arma::mat &opdm);
