/*Copyright (c) 2023 Soubhik Mondal

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


/*! \file BeckeDW1.h
     \brief Calculates 1st derivative of weight parameter w(r_iA) with respect to nucleus A.
 */

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <numgrid.h>
#include "BasisSet.h"
#include <iostream>
#include <fstream>

#ifndef BECKEDW1_H
#define BECKEDW1_H

class BeckeDW1 
{
public:
    size_t num_atoms;
    std::vector<double> center_coordinates_bohrX;
    std::vector<double> center_coordinates_bohrY;
    std::vector<double> center_coordinates_bohrZ;
    std::vector<int> nuc_charges;
    std::vector<Eigen::Vector3d> nuclear_coordinates;
    BeckeDW1(std::vector<Atom>& atoms);

    /**
     * Calculates derivative of weight parameter w(r_iA) with resoect to nucleus A.
     * {r_iA} are grid points.
    */
    void calculate_gridW_derivative(const context_t *context, size_t center_index, 
                                std::map<char, std::map<int, double>>& dB_Weight,  size_t loop_idx);

private:
    /**
     * Calculates the u_ij variable prescribed by Becke in context of Slater-Bragg radii
    */
    double uij( int atom_nuCharge1,  int atom_nuCharge2);
    /**
     * Calculates the a_ij variable which is a function of u_ij
    */
    double aij( int atom_nuCharge1,  int atom_nuCharge2);
    /**
     * Calculates the iterative p(mu) function with Becke hardness of 3.
    */
    double f3(double x);
    /**
     * Calculates auxillary functional value t(mu) for hyberbolic coordinate: mu.
     */
    double tmu(double nu, double tol);
    /**
     * Calculates scalar internuclear distance between nuclei A and B, A and B are indices.
    */
    double R_AB(int A, int B);
    /**
     * Calculates scalar rA-rB, where rX is distance of X th nucleus from the  grid poinr r.
    */
    double rArB( Eigen::Vector3d& r, int A, int B);
    /**
     * Calculates scalar distance between X(i, j, k), Y(i, j, k) coordinates.
    */
    double distanceXY( Eigen::Vector3d& A,  Eigen::Vector3d& B);

    /**
     * Data dumping routines
     * 
    void dumpDataToHDF5(std::string& filename, 
                 std::map<std::pair<int, int>, std::map<char, double>>& dB_PCenter, 
                 std::map<char, std::map<int, double>>& dB_Weight);  
    void dumpDataToCSV(std::string& filename, 
                 std::map<std::pair<int, int>, std::map<char, double>>& dB_PCenter, 
                 std::map<char, std::map<int, double>>& dB_Weight);
   */
};

#endif // BECKEDW1_H
