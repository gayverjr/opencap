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

/*
 * BeckeDW1.cpp
 */
/**
 * References: 
 * 1. Becke, A. D. (1988). A multicenter numerical integration scheme for polyatomic molecules. 
 *    The Journal of chemical physics, 88(4), 2547-2553.
 * 2. Johnson, B. G., Gill, P. M., & Pople, J. A. (1993). The performance of a family of density functional methods. 
 *    The Journal of chemical physics, 98(7), 5612-5626.
*/

#include "BeckeDW1.h"
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <numgrid.h>
#include <bragg.h>
#include "BasisSet.h"
#include <iostream>
#include <fstream>
#include <h5pp/h5pp.h>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


BeckeDW1::BeckeDW1(std::vector<Atom>& atoms):
    num_atoms(atoms.size()),
    center_coordinates_bohrX(num_atoms),
    center_coordinates_bohrY(num_atoms),
    center_coordinates_bohrZ(num_atoms),
    nuc_charges(num_atoms),
    nuclear_coordinates(num_atoms)
 {
    // Set up nuclear charges: Need for bragg radii
    for (size_t i = 0; i < num_atoms; i++) {
            center_coordinates_bohrX[i]=atoms[i].coords[0];
        center_coordinates_bohrY[i]=atoms[i].coords[1];
        center_coordinates_bohrZ[i]=atoms[i].coords[2];
        nuc_charges[i] = atoms[i].Z;
        if (atoms[i].Z == 0)
            nuc_charges[i] = 1; 
    }

    for (int i = 0; i < num_atoms; ++i) {
        nuclear_coordinates[i] = Eigen::Vector3d(center_coordinates_bohrX[i], center_coordinates_bohrY[i], center_coordinates_bohrZ[i]);
    }
}

void BeckeDW1::calculate_gridW_derivative(const context_t *context, size_t center_index, std::map<char, std::map<int, double>>& dB_Weight, size_t loop_idx)
{    
    double tol = 1e-16;
    int num_grid_points = numgrid_get_num_grid_points(context);
    int num_radial_points = numgrid_get_num_radial_grid_points(context);
    int num_angular_grid_points = num_grid_points / num_radial_points;

    double *grid_x_bohr = new double[num_grid_points];
    double *grid_y_bohr = new double[num_grid_points];
    double *grid_z_bohr = new double[num_grid_points];
    double *grid_w = new double[num_grid_points];

    numgrid_get_grid(context, num_atoms, center_index,
                     center_coordinates_bohrX.data(), center_coordinates_bohrY.data(), center_coordinates_bohrZ.data(),
                     nuc_charges.data(), grid_x_bohr, grid_y_bohr, grid_z_bohr,
                     grid_w);
    delete [] grid_w;

    double *radial_grid_r_bohr = new double[num_radial_points];
    double *radial_grid_w = new double[num_radial_points];
    numgrid_get_radial_grid(context,
                            radial_grid_r_bohr,
                            radial_grid_w);
    delete [] radial_grid_r_bohr; 

    double *angular_grid_x_bohr = new double[num_grid_points];
    double *angular_grid_y_bohr = new double[num_grid_points];
    double *angular_grid_z_bohr = new double[num_grid_points];
    double *angular_grid_w = new double[num_grid_points];

    numgrid_get_angular_grid(num_angular_grid_points,
                              angular_grid_x_bohr,
                              angular_grid_y_bohr,
                              angular_grid_z_bohr,
                              angular_grid_w);
    delete [] angular_grid_x_bohr;
    delete [] angular_grid_y_bohr;
    delete [] angular_grid_z_bohr; 

    std::vector<Eigen::Vector3d> coordinates_NC_shifted;
    for (int i = 0; i < num_grid_points; ++i) {
        coordinates_NC_shifted.emplace_back(grid_x_bohr[i], grid_y_bohr[i], grid_z_bohr[i]);
    }

    std::vector<double> wts;
    for (int iRAD = 0; iRAD < num_radial_points; ++iRAD) {
        double wt = 4.0 * M_PI * radial_grid_w[iRAD];
        for (int iANG = 0; iANG < num_angular_grid_points; ++iANG) {
            wts.push_back(wt * angular_grid_w[iANG]);
        }
    }
    delete [] radial_grid_w;
    delete [] angular_grid_w;

    Eigen::MatrixXd P(num_grid_points, num_atoms);
    P.setOnes();
    Eigen::VectorXd Z(num_grid_points); 
    for (int A = 0; A < num_atoms; ++A) {
        for (int B = 0; B < num_atoms; ++B) {
            if (B != A){
                double RAB = R_AB(A, B);
                if (RAB==0.0){
                    std::cout << "Severe warning: Nuclei " << A+1 << " and " << B+1 << " are on top of each other.";
                }
                double aAB = aij(nuc_charges[A], nuc_charges[B]);
                #pragma omp parallel for 
                for (int rCT = 0; rCT < num_grid_points; ++rCT) {
                    double muAB = rArB(coordinates_NC_shifted[rCT], A, B) / RAB;
                    double nuAB = muAB + aAB * (1.0 - muAB * muAB);
                    double smu = 0.50 * (1.0 - f3(f3(f3(nuAB))));
                    P(rCT, A) *= smu;
                    }
                }
            }
        }
    

    Z = P.rowwise().sum();

    std::map<std::pair<int, int>, std::map<char, double>>  dB_PCenter; 
    std::map<std::pair<int, int>, std::map<char, double>>  dB_Z;

    for (int A = 0; A < num_atoms; ++A) {
        for (int i = 0; i < num_grid_points; ++i) {
            dB_Z[std::make_pair(i, A)]['x'] = 0.0;
            dB_Z[std::make_pair(i, A)]['y'] = 0.0;
            dB_Z[std::make_pair(i, A)]['z'] = 0.0;
            dB_PCenter[std::make_pair(i, A)]['x'] = 0.0;
            dB_PCenter[std::make_pair(i, A)]['y'] = 0.0;
            dB_PCenter[std::make_pair(i, A)]['z'] = 0.0;
        }
    }

    for (int B = 0; B < num_atoms; ++B) {
        for (int A = 0; A < num_atoms; ++A) {
            if (A != B) 
            {
                double RAB = R_AB(A, B);
                double aAB = aij(nuc_charges[A], nuc_charges[B]);
                double aBA = aij(nuc_charges[B], nuc_charges[A]);

                #pragma omp parallel for 
                for (int rCT = 0; rCT < num_grid_points; ++rCT) {
                    double x_grid_point = grid_x_bohr[rCT];
                    double y_grid_point = grid_y_bohr[rCT];
                    double z_grid_point = grid_z_bohr[rCT];

                    double muAB = rArB(coordinates_NC_shifted[rCT], A, B) / RAB;
                    double nuAB = muAB + aAB * (1.0 - muAB * muAB);
                    double smuAB = 0.50 * (1.0 - f3(f3(f3(nuAB))));

                    double tmuAB = 0.0;
                    if (std::abs(smuAB) > tol) {
                        tmuAB = (-27.0 / 16.0) * (1.0 - nuAB * nuAB) / smuAB;
                        tmuAB *= (1.0 - f3(f3(nuAB)) * f3(f3(nuAB))) * (1.0 - f3(nuAB) * f3(nuAB));
                    }

                    Eigen::Vector3d uvecB;
                    double factor_distanceB = distanceXY(coordinates_NC_shifted[rCT], nuclear_coordinates[B]);
                    uvecB = Eigen::Vector3d((center_coordinates_bohrX[B]-  x_grid_point)/factor_distanceB ,
                                            (center_coordinates_bohrY[B] - y_grid_point)/factor_distanceB ,
                                            (center_coordinates_bohrZ[B] - z_grid_point)/factor_distanceB);

                    Eigen::Vector3d uvecAB;
                    uvecAB = Eigen::Vector3d((center_coordinates_bohrX[A] - center_coordinates_bohrX[B])/RAB,
                                             (center_coordinates_bohrY[A] - center_coordinates_bohrY[B])/RAB,
                                             (center_coordinates_bohrZ[A] - center_coordinates_bohrZ[B])/RAB);


                    double dBmuABX = (-uvecB.x() + muAB * uvecAB.x()) / RAB;
                    double dBmuABY = (-uvecB.y() + muAB * uvecAB.y()) / RAB;
                    double dBmuABZ = (-uvecB.z() + muAB * uvecAB.z()) / RAB;


                    double dBnuABX = dBmuABX * (1.0 - 2.0 * aAB * muAB);
                    double dBnuABY = dBmuABY * (1.0 - 2.0 * aAB * muAB);
                    double dBnuABZ = dBmuABZ * (1.0 - 2.0 * aAB * muAB);

                    double dBnuBAX = -dBmuABX * (1.0 + 2.0 * aBA * muAB);
                    double dBnuBAY = -dBmuABY * (1.0 + 2.0 * aBA * muAB);
                    double dBnuBAZ = -dBmuABZ * (1.0 + 2.0 * aBA * muAB);

                    /** Eqn. (B8) in BJG 1993 paper. 
                     */

                    if (A==center_index){
                        dB_PCenter[std::make_pair(rCT, B)]['x'] = -tmuAB * P(rCT, A) * dBnuBAX;
                        dB_PCenter[std::make_pair(rCT, B)]['y'] = -tmuAB * P(rCT, A) * dBnuBAY;
                        dB_PCenter[std::make_pair(rCT, B)]['z'] = -tmuAB * P(rCT, A) * dBnuBAZ;
                    }
                    
                    dB_Z[std::make_pair(rCT, B)]['x'] -= tmuAB * P(rCT, A) * dBnuBAX;
                    dB_Z[std::make_pair(rCT, B)]['y'] -= tmuAB * P(rCT, A) * dBnuBAY;
                    dB_Z[std::make_pair(rCT, B)]['z'] -= tmuAB * P(rCT, A) * dBnuBAZ;

                    double tmuBA = 0.0;
                    double nuBA = -muAB + aBA * (1.0 - muAB * muAB);
                    double smuBA = 0.50 * (1.0 - f3(f3(f3(nuBA))));
                    
                    if (std::abs(smuBA) > tol) {
                        tmuBA = (-27.0 / 16.0 * (1.0 - nuBA * nuBA)) / smuBA;
                        tmuBA *= (1.0 - f3(f3(nuBA)) * f3(f3(nuBA))) * (1.0 - f3(nuBA) * f3(nuBA));
                    }
                    dB_Z[std::make_pair(rCT, B)]['x'] += P(rCT, B) * tmuBA * dBnuBAX;
                    dB_Z[std::make_pair(rCT, B)]['y'] += P(rCT, B) * tmuBA * dBnuBAY;
                    dB_Z[std::make_pair(rCT, B)]['z'] += P(rCT, B) * tmuBA * dBnuBAZ;
                }
            }
        }
    }

    for (int i = 0; i < num_grid_points; ++i) {
        dB_Weight['x'][i] = 0.0;
        dB_Weight['y'][i] = 0.0;
        dB_Weight['z'][i] = 0.0;
    }

    for (int B = 0; B < num_atoms; ++B) {
        if (B != center_index) {
            #pragma omp parallel for 
            for (int rCT = 0; rCT < num_grid_points; ++rCT) {
                dB_Weight['x'][rCT] -= (dB_PCenter[std::make_pair(rCT, B)]['x'] / Z[rCT] - P(rCT, center_index) * (dB_Z[std::make_pair(rCT, B)]['x'] / (Z[rCT] * Z[rCT])));
                dB_Weight['y'][rCT] -= (dB_PCenter[std::make_pair(rCT, B)]['y'] / Z[rCT] - P(rCT, center_index) * (dB_Z[std::make_pair(rCT, B)]['y'] / (Z[rCT] * Z[rCT])));
                dB_Weight['z'][rCT] -= (dB_PCenter[std::make_pair(rCT, B)]['z'] / Z[rCT] - P(rCT, center_index) * (dB_Z[std::make_pair(rCT, B)]['z'] / (Z[rCT] * Z[rCT])));
            }
        }
    }

    delete [] grid_x_bohr;
    delete [] grid_y_bohr;
    delete [] grid_z_bohr;

    #pragma omp parallel for 
    for (size_t rCT = 0; rCT < num_grid_points; ++rCT) {
        dB_Weight['x'][rCT] *= wts[rCT];
        dB_Weight['y'][rCT] *= wts[rCT];
        dB_Weight['z'][rCT] *= wts[rCT];
    }
}

double BeckeDW1::uij(int atom_nuCharge1, int atom_nuCharge2) {
    double chi = get_bragg_angstrom(atom_nuCharge1) / get_bragg_angstrom(atom_nuCharge2);
    return (chi - 1) / (chi + 1);
}

double BeckeDW1::aij(int atom_nuCharge1, int atom_nuCharge2) {
    double uij_local = uij(atom_nuCharge1, atom_nuCharge2);
    double value = uij_local / (uij_local * uij_local - 1.0);
    if (value > 0.50) {
        return 0.50;
    } else if (value < -0.50) {
        return -0.50;
    } else {
        return value;
    }
}

double BeckeDW1::f3(double x) {
    return 1.50 * x - 0.5 * x * x * x;
}

double BeckeDW1::tmu(double nu, double tol) {
    double value = 0.0;
    double smu = 0.50 * (1.0 - f3(f3(f3(nu))));
    if (std::abs(smu) > tol){
        value = (-27.0 / 16.0 * (1.0 - nu * nu)) / smu;
        value *= (1.0 - f3(f3(nu)) * f3(f3(nu))) * (1.0 - f3(nu) * f3(nu));}
    return value;
}

double BeckeDW1::R_AB(int A, int B) {
    Eigen::Vector3d CenterA = nuclear_coordinates[A];
    Eigen::Vector3d CenterB = nuclear_coordinates[B];
    return distanceXY(CenterA, CenterB);
}

double BeckeDW1::rArB(Eigen::Vector3d& r, int A, int B) {
    double rA = distanceXY(r, nuclear_coordinates[A]);
    double rB = distanceXY(r, nuclear_coordinates[B]);
    return rA - rB;
}

double BeckeDW1::distanceXY(Eigen::Vector3d& A, Eigen::Vector3d& B) {
    return std::sqrt((A[0] - B[0]) * (A[0] - B[0]) +
                     (A[1] - B[1]) * (A[1] - B[1]) +
                     (A[2] - B[2]) * (A[2] - B[2]));
}
