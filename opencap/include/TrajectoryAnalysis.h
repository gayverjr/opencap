/*Copyright (c) 2020 James Gayvert

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

/*! \file TrajectoryAnalysis.h
     \brief Classes and structs for storing atomic data.
 */

#include <Eigen/Dense>
#include <map>
#include <complex>
#include <cmath>

struct root
{
    double eta;
    std::complex<double> energy;
    Eigen::VectorXcd eigv;
    root(double new_eta, std::complex<double> new_energy, Eigen::VectorXcd new_eigv)
    {eta=new_eta; energy = new_energy; eigv=new_eigv;};
    root()
    {eta=0.0;energy=0.0;};
};


class EigenvalueTrajectory
{
public:
	std::vector<std::complex<double>> uncorrected_energies;
	std::vector<std::complex<double>> corrected_energies;

private:
	root prev;
	std::vector<root> states;

public:
	EigenvalueTrajectory(std::vector<root> initial_states, size_t state_idx);
	void add_state(std::vector<root> new_states);
	void analyze();

};


class CAPHamiltonian
{
public:

private:
	Eigen::MatrixXd CAP_MAT;
	Eigen::MatrixXd ZERO_ORDER_H;
	std::vector<double> eta_list;
	std::vector<EigenvalueTrajectory> trajectories;
	bool do_output;
	size_t nstates;
	std::vector<std::vector<root>> all_roots;
	std::string finame;

public:
	CAPHamiltonian(Eigen::MatrixXd ZERO_ORDER_H,Eigen::MatrixXd CAP_MAT,std::map<std::string, std::string> params,
			std::string filename);
	void run_trajectory();
private:
	void set_eta_list(std::map<std::string, std::string> params);
	void track_states();
	void track_state(size_t state_idx);
	void save_trajectory();
};
