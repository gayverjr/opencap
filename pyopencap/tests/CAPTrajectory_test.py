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
    SOFTWARE.'''

from pyopencap.analysis import CAPHamiltonian
import numpy as np
import os

dest_dir = "../examples/analysis/N2/ref_outputs"
eta_list = np.linspace(0,5000,101) * 1E-5
ref_energy = -109.36195558

def test_eomcc_parse():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"eomcc.out"))
    assert np.isclose(CAPH._H0[0][0],-1.09353e+02)
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"eomcc.out"),irrep="B2g")
    assert np.isclose(CAPH._H0[0][0],-1.09313e+02)
    try:
        CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"eomcc.out"),irrep="B1g")
    except:
        pass

def test_adc_parse():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"adc.out"))
    assert np.isclose(CAPH._W[0][0],-6.3724737e+01)
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"adc.out"),onset="4000")
    assert np.isclose(CAPH._W[0][0],-4.9669293e+01)
    try:
        CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"adc.out"),onset="1000")
    except:
        pass

def test_opencap_parse():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"n2_opencap.out"))
    assert np.isclose(CAPH._W[0][0],-4.28708660e+01)

def test_eta_opt():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"n2_opencap.out"))
    CAPH.run_trajectory(eta_list)
    traj = CAPH.track_state(1,tracking="overlap")
    uc_energy_au,uc_eta = traj.find_eta_opt(units="au")
    corr_energy_au,corr_eta = traj.find_eta_opt(corrected=True,units="au")
    uc_energy_ev,uc_eta = traj.find_eta_opt(ref_energy=ref_energy,units="eV")
    corr_energy_ev,corr_eta = traj.find_eta_opt(ref_energy=ref_energy,corrected=True,units="eV")
    assert np.isclose(np.real(uc_energy_au),-109.265006)
    assert np.isclose(np.real(corr_energy_au), -109.26652)
    assert np.isclose(np.real(uc_energy_ev),2.6381)
    assert np.isclose(np.real(corr_energy_ev),2.59686)
    traj = CAPH.track_state(1,tracking="energy")
    uc_energy_au = traj.get_energy(uc_eta)
    corr_energy_au = traj.get_energy(corr_eta,corrected=True)
    uc_energy_ev = traj.get_energy(uc_eta,ref_energy=ref_energy,units="eV")
    corr_energy_ev = traj.get_energy(corr_eta,ref_energy=ref_energy,corrected=True,units="eV")
    assert np.isclose(np.real(uc_energy_au),-109.265006)
    assert np.isclose(np.real(corr_energy_au), -109.26652)
    assert np.isclose(np.real(uc_energy_ev),2.6381)
    assert np.isclose(np.real(corr_energy_ev),2.59686)

def test_exclude_states():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"n2_opencap.out"))
    CAPH.run_trajectory(eta_list,exclude_states=[3])
    assert CAPH.H0.shape == CAPH.W.shape == (4,4)
    assert np.isclose(CAPH.H0[3][3], -1.089933e+02)
    assert np.isclose(CAPH.W[3][3], -1.65280e+02)

def test_include_states():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"n2_opencap.out"))
    CAPH.run_trajectory(eta_list,include_states=[0,1,2,4])
    assert CAPH.H0.shape == CAPH.W.shape == (4,4)
    assert np.isclose(CAPH.H0[3][3], -1.089933e+02)
    assert np.isclose(CAPH.W[3][3], -1.65280e+02)

def test_export():
    CAPH = CAPHamiltonian(output=os.path.join(dest_dir,"n2_opencap.out"))
    CAPH.export("test.out")
    CAPH2 = CAPHamiltonian(output="test.out")
    assert np.isclose(CAPH._H0[0][0],CAPH2._H0[0][0])
    os.remove("test.out")
