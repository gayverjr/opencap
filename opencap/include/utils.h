 /*! \file utils.h
     \brief Helper functions used by various pieces of the software.
 */
#include "Shell.h"
#include <Eigen/Dense>
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
void fill_LT(std::vector<double> matrix_elements, Eigen::MatrixXd &opdm);
/*! Populates matrix given full list of elements in row major order
 * \param matrix_elements: list of all matrix elements in row-major order
 * \param opdm: matrix to be populated
*/
void fill_mat(std::vector<double> &matrix_elements, Eigen::MatrixXd &opdm);
