import os
import numpy as np
import functools
from numpy import linalg as LA
import matplotlib.pyplot as plt

ccsd_energy = -109.35451902
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
    

irrep="B2g"
with open('caspt2.out', 'r') as file :
    filedata = file.readlines()

for i in range(0,len(filedata)):
	line = filedata[i]
	if "Starting Perturbative CAP calculation." in line:
		break
for j in range(i,len(filedata)):
	line=filedata[j]
	if "Solving for CAP-EOM-CCSD " + irrep + " transitions" in line:
		break
num_roots=int(filedata[j+1].split()[4])
start=j+5
cap_mat=[]
for i in range(start,start+num_roots):
	l1 = filedata[i].split()
	l1= [float(x) for x in l1]
	cap_mat+=l1
cap_mat = np.reshape(cap_mat,(num_roots,num_roots))


for q in range(0,len(cap_mat)):
    for w in range(q,len(cap_mat)):
        cap_mat[w][q]=cap_mat[q][w]

start=i+2
H_0=[]
for i in range(start,start+num_roots):
	l1 = filedata[i].split()
	l1= [float(x) for x in l1]
	H_0+=l1
H_0 = np.reshape(H_0,(num_roots,num_roots))	


guess=2.5
eta_list = np.linspace(0,2000,201)
eta_list = eta_list * 1E-5
all_roots=[]
for i in range(0,len(eta_list)):
    eta=eta_list[i]
    roots=[]
    fullH = H_0 +1.0j * eta * cap_mat
    eigv,eigvc=LA.eig(fullH)
    for eig in eigv:
        E = (eig - ccsd_energy) * au2eV
        roots.append(root(E,eta))
        all_roots.append(root(E,eta))
    if i==0:
        traj=trajectory(roots,guess)
    else:
        traj.add_state(roots)


re_traj = []
im_traj = []
print("Printing energies")
energies=[]
for root in all_roots:
    re_traj.append(np.real(root.energy))
    print(root.energy)
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
#plt.plot(corr_re,corr_im,'-bo')
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



