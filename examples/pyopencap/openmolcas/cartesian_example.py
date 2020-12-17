import pyopencap
import numpy as np
from pandas import DataFrame
import h5py
from pyscf import gto, scf, ci, ao2mo

#Change these lines to suit your system
##########################################
ref_energy = -109.35814085
au2eV= 27.2113961
guess = 2.2
eta_list = np.linspace(0,500,51)
eta_list = eta_list * 1E-5
RASSI_FILE = "../../opencap/cart.rassi.h5"
OUTPUT_FILE = "../../opencap/cart.out"
nbasis=140
##########################################

sys_dict = {"molecule": "molcas_rassi",
"basis_file": RASSI_FILE}

cap_dict = {
    "cap_type": "box",
        "cap_x":"2.76",
            "cap_y":"2.76",
            "cap_z":"4.88",
            "Radial_precision": "14",
            "angular_points": "110"
}

es_dict = {"method" : "ms-caspt2",
    "molcas_output":OUTPUT_FILE ,
        "rassi_h5":RASSI_FILE,
}

f = h5py.File(RASSI_FILE, 'r')
dms = f["SFS_TRANSITION_DENSITIES"]
overlap_mat = np.array(f["AO_OVERLAP_MATRIX"])
overlap_mat= np.reshape(overlap_mat,(nbasis,nbasis))

# Method 1: Use renormalize function after reading in data
s = pyopencap.System(sys_dict)
s.check_overlap_mat(overlap_mat,"openmolcas",RASSI_FILE)
pc = pyopencap.CAP(s,cap_dict,10,"openmolcas")
pc.read_data(es_dict)
pc.compute_ao_cap()
pc.renormalize()
pc.compute_perturb_cap()
mat= pc.get_perturb_cap()
h0 = pc.get_H()


# Method 2: Use renormalize_cap to renormalize after reading in DMs separately
pc = pyopencap.CAP(s,cap_dict,10,"openmolcas")
for i in range(0,10):
    for j in range(i,10):
        dm1 = 0.5*np.reshape(dms[i][j],(nbasis,nbasis))
        pc.add_tdms(dm1,dm1,i,j,"openmolcas",RASSI_FILE)
        if i!=j:
            pc.add_tdms(dm1,dm1,j,i,"openmolcas",RASSI_FILE)
pc.compute_ao_cap()
pc.renormalize_cap(overlap_mat,"openmolcas",RASSI_FILE)
pc.compute_perturb_cap()
mat=pc.get_perturb_cap()



### now time for trajectories
import os
import functools
from numpy import linalg as LA
import matplotlib.pyplot as plt


# a root is a single eigenvalue of the cap hamiltonian at a particular value of eta
@functools.total_ordering
class root():
    def __init__(self, energy, eta):
        self.eta = eta
        self.energy = energy
    
    def __lt__(self, other):
        return self.eta < other.eta
    
    def __eq__(self, other):
        return self.eta == other.eta and self.eta == other.eta

# a trajectory is a group of eigenvalues over a range of eta values, grouped by proximity to an initial guess
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
    
    # picks out the state from the list of states whose energy is closest to the previous entry in the trajectory
    def add_state(self,states):
        min=500
        cur=-1
        for st in states:
            if np.absolute(st.energy-self.last.energy)<min:
                cur=st
                min=np.absolute(st.energy-self.last.energy)
        self.last = cur
        self.states.append(cur)
    # applies first order correciton
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
# A previous run through of this script showed the resonance trajectory starting near 2.2eV, so that'll be our initial guess
all_roots=[]
# diagonalize over range of eta values and generate trajectories
for i in range(0,len(eta_list)):
    eta=eta_list[i]
    roots=[]
    fullH = H_0 +1.0j * eta * cap_mat
    eigv,eigvc=LA.eig(fullH)
    for eig in eigv:
        E = (eig - ref_energy) * au2eV
        roots.append(root(E,eta))
        all_roots.append(root(E,eta))
    if i==0:
        traj=trajectory(roots,guess)
    else:
        traj.add_state(roots)

# first lets plot everything
re_traj = []
im_traj = []
energies=[]
for root in all_roots:
    re_traj.append(np.real(root.energy))
    im_traj.append(np.imag(root.energy))
    energies.append(root.energy)
plt.title("Eigenvalue trajectories")
plt.xlabel("Re(E)[eV]")
plt.ylabel("Im(E)[eV]")
plt.plot(re_traj,im_traj,'ro')
plt.show()

# lets get the corrected trajectory
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
# plot uncorrected and corrected trajectory
plt.title("Resonance trajectory")
plt.plot(re_traj,im_traj,'-ro',label="Uncorrected trajectory")
plt.plot(corr_re,corr_im,'-bo',label="Corrected trajectory")
plt.xlabel("Re(E)[eV]")
plt.ylabel("Im(E)[eV]")
plt.legend()
plt.show()

# plot derivative, find stationary point on uncorrected trajectory
derivs=list(np.array(eta_list)*np.absolute(np.gradient(uc_energies)/np.gradient(eta_list)))
plt.plot(eta_list,derivs)
plt.title("Uncorrected derivative")
plt.show()
sorted_derivs = sorted(derivs[5:])
points = []
etas = []
for i in range(0,5):
    points.append(uc_energies[derivs.index(sorted_derivs[i])])
    etas.append(eta_list[derivs.index(sorted_derivs[i])])
print("Uncorrected:")
print(points)
print(sorted_derivs[:5])
print(etas)

# plot derivative, find stationary point on corrected trajectory
derivs=list(np.array(eta_list)*np.absolute(np.gradient(corr_energies)/np.gradient(eta_list)))
plt.plot(eta_list,derivs)
plt.title("Corrected derivative")
plt.show()
sorted_derivs = sorted(derivs[5:])
points = []
etas = []
for i in range(0,5):
    points.append(corr_energies[derivs.index(sorted_derivs[i])])
    etas.append(eta_list[derivs.index(sorted_derivs[i])])
print("Corrected:")
print(points)
print(sorted_derivs[:5])
print(etas)





