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

#include "TrajectoryAnalysis.h"

#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>


#include "utils.h"

void CAPHamiltonian::set_eta_list(std::map<std::string, std::string> params)
{
	double eta_step;
	size_t nsteps;
	std::stringstream stepss(params["eta_step"]);
	std::stringstream nstepss(params["nsteps"]);
	std::istringstream save_traj(params["save_trajectory"]);
	stepss >> eta_step;
	nstepss >> nsteps;
	save_traj >> std::boolalpha >> do_output;
	eta_step*=1E-5;
	for(size_t i=0;i<=nsteps;i++)
		eta_list.push_back(i*eta_step);
}

CAPHamiltonian::CAPHamiltonian(Eigen::MatrixXd h0,Eigen::MatrixXd cap,std::map<std::string,
		std::string> params,std::string filename)
{
	CAP_MAT = cap;
	ZERO_ORDER_H = h0;
	nstates = h0.cols();
	set_eta_list(params);
	finame = filename;
}

void CAPHamiltonian::track_state(size_t state_idx)
{
	std::cout << "Results for state:" << state_idx << std::endl;
	std::cout << "----------------------" << std::endl;
	EigenvalueTrajectory traj(all_roots[0],state_idx);
	for(size_t i=1;i<all_roots.size();i++)
		traj.add_state(all_roots[i]);
	traj.analyze();
	std::cout << "----------------------" << std::endl;
	trajectories.push_back(traj);
}

void CAPHamiltonian::track_states()
{
	std::cout << "---------------------" << std::endl;
	std::cout << "Trajectory Analysis:" << std::endl;
	std::cout << "--------------------" << std::endl;
	for(size_t i=0;i<nstates;i++)
		track_state(i);
}

void CAPHamiltonian::save_trajectory()
{
	std::ofstream trajFile;
	std::string data_finame = finame + ".data";
	trajFile.open(data_finame);
    trajFile << std::left << std::setw(18) << std::setfill(' ')
    << "State" << std::left << std::setw(18) << std::setfill(' ')
    << "eta" << std::left << std::setw(18) << std::setfill(' ')
    << "Uncorr_Re" << std::left << std::setw(18) << std::setfill(' ')
    << "Uncorr_Im" << std::left << std::setw(18) << std::setfill(' ');
    trajFile << "Corr_Re" << std::left << std::setw(18) << std::setfill(' ') << "Corr_Im" << std::endl;
    trajFile.close();
}

void CAPHamiltonian::run_trajectory()
{
	size_t nstates = CAP_MAT.rows();
	for(size_t i=0;i<eta_list.size();i++)
	{
		std::vector<root> roots;
		Eigen::MatrixXcd CAPH(nstates,nstates);
		CAPH.real() = ZERO_ORDER_H;
		CAPH.imag() = eta_list[i]*CAP_MAT;
		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
		ces.compute(CAPH);
		for(size_t j=0;j<nstates;j++)
		{
			std::complex<double> energy =  ces.eigenvalues()[j];
			Eigen::VectorXcd vec = ces.eigenvectors().col(j);
			root r1(eta_list[i],energy,vec);
			roots.push_back(r1);
		}
		all_roots.push_back(roots);
	}
	track_states();
	if(do_output)
		save_trajectory();
}

EigenvalueTrajectory::EigenvalueTrajectory(std::vector<root> initial_states, size_t state_idx)
{
	prev = initial_states[state_idx];
	states.push_back(prev);
	uncorrected_energies.push_back(prev.energy);
	corrected_energies.push_back(prev.energy);
}

void EigenvalueTrajectory::add_state(std::vector<root> new_states)
{
	double max_overlap=0.0;
	root best=new_states[0];
	for(auto cur: new_states)
	{
		double overlap = abs(prev.eigv.dot(cur.eigv));
		if(overlap>max_overlap)
		{
			best = cur;
			max_overlap = overlap;
		}
	}
	states.push_back(best);
	prev = best;
	uncorrected_energies.push_back(best.energy);
}

void EigenvalueTrajectory::analyze()
{
	// get corrections, find minimum of logarithmic velocity of uncorrected trajectory

	size_t num_points = uncorrected_energies.size();
	std::vector<double> grad;
	std::vector<double> corr_grad;
	double step = states[1].eta;

	grad.push_back(0.0);
	corrected_energies.push_back(uncorrected_energies[0]);
	for (int i=1;i<uncorrected_energies.size()-1;i++)
	{
		std::complex<double> log_velo = states[i].eta*(uncorrected_energies[i+1]-uncorrected_energies[i-1])/(2.0*step);
	    grad.push_back(std::abs(log_velo));
	    corrected_energies.push_back(uncorrected_energies[i]-log_velo);
	}
	// last point
	std::complex<double> last_deriv = states[num_points-1].eta*(uncorrected_energies[num_points-1]-uncorrected_energies[num_points-2])/step;
	grad.push_back(std::abs(last_deriv));
	corrected_energies.push_back(uncorrected_energies[num_points-1]-last_deriv);

	corr_grad.push_back(0.0);
	for (int i=1;i<corrected_energies.size();i++)
	    corr_grad.push_back(states[i].eta*std::abs((corrected_energies[i+1]-corrected_energies[i-1])/(2.0*step)));

	double min_deriv=1000;
	size_t min_idx = 0;
	double min_deriv_corr=1000;
	size_t min_idx_corr;
	for(size_t i=10;i<grad.size()-1;i++)
	{
		if(grad[i]<min_deriv)
		{
			min_deriv = grad[i];
			min_idx = i;
		}
		if(corr_grad[i]<min_deriv_corr)
		{
			min_deriv_corr=corr_grad[i];
			min_idx_corr = i;
		}
	}


	//uncorrected results
	std::cout << "Results from uncorrected trajectory:" << std::endl;
	std::cout << "Uncorrected energy: " << std::noshowpos << uncorrected_energies[min_idx] << std::endl;
	std::cout << "Eta: " << states[min_idx].eta << " Logarithmic velocity: " << min_deriv << std::endl;
	uc_opt_eta = states[min_idx].eta;
	uc_opt = uncorrected_energies[min_idx];

	// corrected results
	std::cout << "Results from corrected trajectory:" << std::endl;
	std::cout << "Corrected energy: " << std::noshowpos << corrected_energies[min_idx_corr] << std::endl;
	std::cout << "Eta: " << states[min_idx_corr].eta << " Logarithmic velocity: " << min_deriv_corr << std::endl;
	corr_opt_eta = states[min_idx].eta;
	corr_opt = uncorrected_energies[min_idx];
}
