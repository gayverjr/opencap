'''Copyright (c) 2022 James Gayvert, Soubhik Mondal
    
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
''' 
MR-CISD+pople correction on N2- with aug-cc-pvtz+3s3p3d basis set.
'''


import pyopencap
from pyopencap.analysis import CAPHamiltonian,colparser
import numpy as np
import matplotlib.pyplot as plt



nstates = 4
molden_dict = {"basis_file":"data_files/molden_mo_mc.sp",
"molecule": "molden"}
s = pyopencap.System(molden_dict)
cap_dict = {"cap_type": "voronoi","r_cut": "3.00"}

pc = pyopencap.CAP(s,cap_dict,nstates)
# Parsing of Columbus data is done through pyopencap.analysis.colparser object
parser = colparser('data_files/molden_mo_mc.sp', 'data_files/tranls')
print(parser.mo_summary())
H0 = parser.get_H0(filename='data_files/ciudgsm')
for i in range(0,nstates):
    for j in range(i,nstates):
        if i==j:
            # Indices start from 0 in pyopencap, but from 1 in Columbus file names
            dm1_ao = parser.sdm_ao(i+1,data_dir='data_files')
            pc.add_tdm(dm1_ao,i,j,'molden')
        else:
            # Indices start from 0 in pyopencap, but from 1 in Columbus file names
            dm1_ao = parser.tdm_ao(i+1, j+1,data_dir='data_files')
            pc.add_tdm(dm1_ao,i,j,'molden')
            pc.add_tdm(dm1_ao.conj().T,j,i,'molden')
pc.compute_projected_cap()
W=pc.get_projected_cap()

CAPH = CAPHamiltonian(H0=H0,W=W)
# eci+pople energy from data_files/ciudgsm_neutral
ref_energy = -109.379094090794
eta_list = np.linspace(0,2000,1001) * 1E-5
CAPH.run_trajectory(eta_list)
plt.plot(np.real(CAPH.energies_ev(ref_energy=ref_energy)),np.imag(CAPH.energies_ev(ref_energy=ref_energy)),'ro',label='Uncorrected Trajectory')
plt.show()

traj = CAPH.track_state(1, tracking="overlap")
uc_energies = traj.energies_ev(ref_energy=ref_energy)
corr_energies = traj.energies_ev(ref_energy=ref_energy,corrected=True)
uc_energy, eta_opt = traj.find_eta_opt(start_idx=20)
uc_energy = (uc_energy-ref_energy)*27.2114
corr_energy, corr_eta_opt = traj.find_eta_opt(corrected=True,start_idx=20)
corr_energy = (corr_energy-ref_energy)*27.2114
print("Uncorrected:")
print(uc_energy)
print(eta_opt)
print("Corrected:")
print(corr_energy)
print(corr_eta_opt)
plt.plot(np.real(uc_energies),np.imag(uc_energies),'-mo',label="Uncorrected")
plt.plot(np.real(corr_energies),np.imag(corr_energies),'-bo',label="Corrected")
plt.plot(np.real(uc_energy),np.imag(uc_energy),'g*',markersize=20)
plt.plot(np.real(corr_energy),np.imag(corr_energy),'r*',markersize=20)
plt.legend()
plt.show()


