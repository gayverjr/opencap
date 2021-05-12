import os
import numpy as np
import functools
from numpy import linalg as LA
import matplotlib.pyplot as plt
import argparse

#########################################
# Alter these values to suit your purposes
ref_energy = -109.36009153
guess = 2.2
eta_list = np.linspace(0,5000,101)
###########################
au2eV= 27.2113961
eta_list = eta_list * 1E-5

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


# read in data from opencap output file
parser = argparse.ArgumentParser()
parser.add_argument('pos_arg', type=str,
                    help='Name of OpenCAP output file')
args = parser.parse_args()
with open(args.pos_arg, 'r') as file :
    filedata = file.readlines()

idx=-1
for i in range(0,len(filedata)):
    line = filedata[i]
    if "Printing out matrices required for Projected CAP calculation." in line:
        idx=i
num_roots=int(filedata[idx+1].split()[-1])
start=idx+3
# zeroth order hamiltonian first
H_0=[]
for i in range(start,start+num_roots):
    l1 = filedata[i].split()
    l1= [float(x) for x in l1]
    H_0+=l1
H_0 = np.reshape(H_0,(num_roots,num_roots))

start2=start+num_roots+1
# now the cap matrix
cap_mat=[]
for i in range(start2,start2+num_roots):
    l1 = filedata[i].split()
    l1= [float(x) for x in l1]
    cap_mat+=l1
cap_mat= np.reshape(cap_mat,(num_roots,num_roots))

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

