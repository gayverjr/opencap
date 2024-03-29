'''Copyright (c) 2021 James Gayvert
    
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
    '''
    
# Analyze OpenCAP output file, view eigenvalue trajectories
from pyopencap.analysis import CAPHamiltonian
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Process OpenCAP output file.')
parser.add_argument('output_file',help='OpenCAP output file')
args = parser.parse_args()


CAPH = CAPHamiltonian(output=args.output_file)
eta_list = np.linspace(0,5000,101)
eta_list = np.around(eta_list * 1E-5,decimals=5)
CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.total_energies),np.imag(CAPH.total_energies),'ro')
plt.show()

for i in range(0,CAPH.nstates):
	traj = CAPH.track_state(i,tracking="overlap")
	uc_energy, eta_opt = traj.find_eta_opt(start_idx=10)
	corr_energy, corr_eta_opt = traj.find_eta_opt(corrected=True,start_idx=10)
	plt.plot(np.real(traj.uncorrected_energies),np.imag(traj.uncorrected_energies),'-ro')
	plt.plot(np.real(traj.corrected_energies),np.imag(traj.corrected_energies),'-bo')
	plt.plot(np.real(corr_energy),np.imag(corr_energy),"g*",markersize=20)
	plt.plot(np.real(uc_energy),np.imag(uc_energy),"g*",markersize=20)
	plt.show()

# Now in eV with reference energy
ref_energy = 0.0 
for i in range(0,CAPH.nstates):  
    plt.title(i)  
    traj = CAPH.track_state(i,tracking="energy",correction="density")
    uc_energies = traj.energies_ev(ref_energy=ref_energy)
    corr_energies = traj.energies_ev(ref_energy=ref_energy,corrected=True)
    uc_energy, eta_opt = traj.find_eta_opt(start_idx=5,end_idx=-20)
    uc_energy = (uc_energy-ref_energy)*27.2114
    corr_energy, corr_eta_opt = traj.find_eta_opt(corrected=True,start_idx=10)
    corr_energy = (corr_energy-ref_energy)*27.2114
    plt.plot(np.real(uc_energies),np.imag(uc_energies),'-ro',label="Uncorrected")
    plt.plot(np.real(corr_energies),np.imag(corr_energies),'-bo',label="Corrected")
    plt.plot(np.real(uc_energy),np.imag(uc_energy),'g*',markersize=20)
    plt.plot(np.real(corr_energy),np.imag(corr_energy),'g*',markersize=20)
    plt.legend()
    plt.show()