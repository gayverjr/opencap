import pyopencap
import numpy as np
from pandas import DataFrame
import h5py
from pyscf import gto, scf, ci, ao2mo

sys_dict = {"geometry":    '''N  0  0   1.039
                             N  0  0   -1.039
                            X 0  0   0.0''',
            "molecule" : "inline",
            "basis_file":"molcas_bas.bas"}

cap_dict = {
            "cap_type": "box",
            "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

es_dict = {"method" : "ms-caspt2",
           "molcas_output":"anion_xms.out",
           "rassi_h5":"anion_xms.rassi.h5",
}

f = h5py.File('anion_xms.rassi.h5', 'r')
arr = f["SFS_TRANSITION_DENSITIES"]
arr2 = np.array(f["AO_OVERLAP_MATRIX"])
arr2 = np.reshape(arr2,(119,119))

#read data
s = pyopencap.System(sys_dict)
s.check_overlap_mat(arr2,"openmolcas","anion_xms.rassi.h5")
pc = pyopencap.Projected_CAP(s,cap_dict,10,"openmolcas")
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.compute_projected_cap()
mat=pc.get_projected_cap()
h0 = pc.get_H()
#print(DataFrame(h0).to_string(index=False, header=False))
#print(DataFrame(mat).to_string(index=False, header=False))


# separate alpha beta
pc = pyopencap.Projected_CAP(s,cap_dict,10,"openmolcas")
for i in range(0,10):
    for j in range(i,10):
        arr1 = 0.5*np.reshape(arr[i][j],(119,119))
        pc.add_tdms(arr1,arr1,i,j,"openmolcas","anion_xms.rassi.h5")
        if i!=j:
            pc.add_tdms(arr1,arr1,j,i,"openmolcas","anion_xms.rassi.h5")
pc.compute_ao_cap()
pc.compute_projected_cap()
mat=pc.get_projected_cap()
#print(DataFrame(mat).to_string(index=False, header=False))


# spin traced
pc = pyopencap.Projected_CAP(s,cap_dict,10,"openmolcas")
for i in range(0,10):
    for j in range(i,10):
        arr1 = np.reshape(arr[i][j],(119,119))
        pc.add_tdm(arr1,i,j,"openmolcas","anion_xms.rassi.h5")
        if i!=j:
            pc.add_tdm(arr1,j,i,"openmolcas","anion_xms.rassi.h5")
pc.compute_ao_cap()
pc.compute_projected_cap()
mat=pc.get_projected_cap()
print(DataFrame(mat).to_string(index=False, header=False))


### now time for trajectories
import os
import functools
from numpy import linalg as LA
import matplotlib.pyplot as plt
xms_caspt2_energy = -109.35714881
E_0 = xms_caspt2_energy
au2eV= 27.2113961
@functools.total_ordering
class root():
    def __init__(self, energy, eta):
        self.eta = eta
        self.energy = energy
    
    def __lt__(self, other):
        return self.eta < other.eta
    
    def __eq__(self, other):
        return self.eta == other.eta and self.eta == other.eta

class trajectory():
    def __init__(self,states,guess):
        min=500
        cur=-1
        for st in states:
            if np.absolute(st.energy-guess)<min:
                cur=st
                min=np.absolute(st.energy-guess)
        self.last=cur
        self.states=[cur]
    
    def add_state(self,states):
        min=500
        cur=-1
        for st in states:
            if np.absolute(st.energy-guess)<min:
                cur=st
                min=np.absolute(st.energy-self.last.energy)
        self.states.append(cur)
    
    def get_corrections(self):
        energies=[]
        etas=[]
        for st in self.states:
            energies.append(st.energy)
            etas.append(st.eta)
        derivs=list(np.gradient(energies)/np.gradient(etas))
        for i in range(0,len(self.states)):
            self.states[i].corr_energy=self.states[i].energy-derivs[i]*self.states[i].eta

H_0 = h0
cap_mat = mat
guess=3
eta_list = np.linspace(0,2000,201)
eta_list = eta_list * 1E-5
all_roots=[]
for i in range(0,len(eta_list)):
    eta=eta_list[i]
    roots=[]
    fullH = H_0 +1.0j * eta * cap_mat
    eigv,eigvc=LA.eig(fullH)
    for eig in eigv:
        E = (eig - E_0) * au2eV
        roots.append(root(E,eta))
        all_roots.append(root(E,eta))
    if i==0:
        traj=trajectory(roots,guess)
    else:
        traj.add_state(roots)
re_traj = []
im_traj = []
energies=[]
for root in all_roots:
    re_traj.append(np.real(root.energy))
    im_traj.append(np.imag(root.energy))
    energies.append(root.energy)
plt.plot(re_traj,im_traj,'ro')
plt.show()


traj.get_corrections()
re_traj = []
im_traj = []
corr_re=[]
corr_im=[]
uc_energies=[]
corr_energies=[]
for root in traj.states:
    uc_energies.append(root.energy)
    re_traj.append(np.real(root.energy))
    im_traj.append(np.imag(root.energy))
    corr_re.append(np.real(root.corr_energy))
    corr_im.append(np.imag(root.corr_energy))
    corr_energies.append(root.corr_energy)
plt.plot(re_traj,im_traj,'-ro')
plt.plot(corr_re,corr_im,'-bo')
plt.show()

derivs=list(np.absolute(np.gradient(uc_energies)/np.gradient(eta_list)))
plt.plot(eta_list,derivs)
plt.title("UC deriv")
plt.show()
sorted_derivs = sorted(derivs)
points = []
etas = []
for i in range(0,5):
    points.append(uc_energies[derivs.index(sorted_derivs[i])])
    etas.append(eta_list[derivs.index(sorted_derivs[i])])
print("Uncorrected:")
print(points)
print(sorted_derivs[:5])
print(etas)


derivs=list(np.absolute(np.gradient(corr_energies)/np.gradient(eta_list)))
plt.plot(eta_list,derivs)
plt.title("Corr deriv")
plt.show()
sorted_derivs = sorted(derivs)
points = []
etas = []
for i in range(0,5):
    points.append(corr_energies[derivs.index(sorted_derivs[i])])
    etas.append(eta_list[derivs.index(sorted_derivs[i])])
print("Corrected together:")
print(points)
print(sorted_derivs[:5])
print(etas)


plt.plot(eta_list,corr_re)
plt.title("Real part")
plt.show()
derivs=list(np.absolute(np.gradient(corr_re)/np.gradient(eta_list)))
plt.plot(eta_list,derivs)
plt.title("Real deriv")
plt.show()
sorted_derivs = sorted(derivs)
points = []
etas = []
for i in range(0,5):
    points.append(corr_re[derivs.index(sorted_derivs[i])])
    etas.append(eta_list[derivs.index(sorted_derivs[i])])
print("Real part:")
print(points)
print(sorted_derivs[:5])
print(etas)


plt.plot(eta_list,corr_im)
plt.title("Imag part")
plt.show()
derivs=list(np.absolute(np.gradient(corr_im)/np.gradient(eta_list)))
plt.plot(eta_list,derivs)
plt.title("Imag deriv")
plt.show()
sorted_derivs = sorted(derivs)
points = []
etas = []
for i in range(0,5):
    points.append(corr_im[derivs.index(sorted_derivs[i])])
    etas.append(eta_list[derivs.index(sorted_derivs[i])])
print("Imaginary part:")
print(points)
print(sorted_derivs[:5])
print(etas)



