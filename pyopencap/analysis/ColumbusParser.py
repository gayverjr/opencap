import os, sys
import re
import numpy as np
import pandas as pd
pd.set_option("display.precision", 16)
from pathlib2 import Path
from scipy.linalg import block_diag
from pandas import DataFrame as df
import pyopencap
from pyopencap.analysis import CAPHamiltonian



class MO:
    def __init__(self,line,total_idx):
        irrep_label = line.split('=')[1].lower()
        irrep_label = ''.join(c for c in irrep_label if c.isalnum())
        first_letter = irrep_label.find(next(filter(str.isalpha, irrep_label)))
        self.irrep = irrep_label[first_letter:]
        self.orb_num_symm = irrep_label[:first_letter]
        self.orb_num_total = total_idx
        self.name = str(self.orb_num_symm) + self.irrep
        self.coeffs = []
    def add_coeff(self,coeff):
        self.coeffs.append(coeff)


#Columbus parser class
class colparser:

    def set_mo_coeff(self,ordering=None):
        if ordering is not None:
            raise NotImplementedError("Custom ordering NYI.")
        else:
            mo_coeff = np.zeros((self.nbft,self.nbft))
            for i, mo in enumerate(self.mos):
                mo_coeff[:,i] = mo.coeffs
            self.mo_coeff = mo_coeff


    def read_molden(self,molden_file):
        _SEC_REGEX = re.compile(r'\[[^]]+\]')
        irrep_dict = {}
        with open(molden_file,'r') as f:
            line = f.readline()
            while '[MO]' not in line:
                line = f.readline()
            cur_mo = None
            mos = []
            line = f.readline()
            new_sec = _SEC_REGEX.match(line)
            while new_sec is None and len(line.strip())!=0:
                if 'sym=' in line.lower():
                    if cur_mo is not None:
                        mos.append(cur_mo)
                    cur_mo = MO(line,len(mos)+1)
                    if cur_mo.irrep in irrep_dict:
                        irrep_dict[cur_mo.irrep] = irrep_dict[cur_mo.irrep]+1
                    else:
                        irrep_dict[cur_mo.irrep] = 1
                elif 'ene=' in line.lower():
                    pass
                elif 'spin' in line.lower():
                    pass
                elif 'occup' in line.lower():
                    pass
                elif len(line.strip())!=0:
                    coeff = float(line.split()[1])
                    cur_mo.add_coeff(coeff)
                line = f.readline()
                new_sec = _SEC_REGEX.match(line)
            mos.append(cur_mo)
        return mos


    def __init__(self, molden_file, tranls, nstates):
        
        self.nstates = nstates

        rfile = open(tranls, 'r')
        while True:
            try:
                line=next(rfile)
            except StopIteration:
                break
            if 'transformation information:' in line:
                line=next(rfile)
                line=next(rfile)
                self.nbft = int(re.findall(r'\d+', line)[0])
                continue

            if 'symmetry blocking information:' in line:
                line=next(rfile)
                line=next(rfile)
                self.slabel = ([str(x) for x in line.split()[2:len(line)]])
                line=next(rfile)
                self.nbpsy = ([int(x) for x in line.split()[2:len(line)]])
                line=next(rfile)
                self.nmpsy = ([int(x) for x in line.split()[2:len(line)]])
                line=next(rfile)
                self.nfcpsy = ([int(x) for x in line.split()[2:len(line)]])
                break

        self.mos = self.read_molden(molden_file)
        self.set_mo_coeff()


    def irrep_info(self):
        pformat = r' Total number of basis functions: {}' \
        ' \n Symmetry Label: {}' \
        ' \n NBPSY: {}'' \n NMPSY: {}' ' \n NFCPSY: {} \n'
        return pformat.format(self.nbft, self.slabel, self.nbpsy, self.nmpsy, self.nfcpsy)

    
    #Gives DM in MO basis
    def dm_from_iwfmt(self, finame, state_dm=False):
        if state_dm:occ = 2.0
        else: occ = 0.0


        with open(finame,'r') as f:
            lines = f.readlines()

        dm1_mo = np.zeros((self.nbft,self.nbft))
        in_data = False
        sym_store = []
        symm_indices = []
        asymm_store = {}
        idx = 4

        while True:
            try:
                line = lines[idx]
            except IndexError:
                break
            if ' 0.000000000000E+00' in line:
                if len(line.split()) == 1:
                    words = last_line.split()
                    (num,lab1,ibvtyp,itypea,itypeb,ifmt,last,nipv) = [int(word) for word in words]
                    if itypea==0 and itypeb==7:
                        in_data=True
                        sym = 1
                    elif itypea==2 and itypeb==9:
                        in_data=True
                        sym = -1
            elif '    ' in line:
                pass
            elif in_data:
                words = line.split()
                val = float(words[0])
                i = int(words[1])-1
                j = int(words[2])-1
                if sym==1:
                    sym_store.append(val)
                    symm_indices.append((i,j))
                elif sym==-1:
                    asymm_store[(i,j)] = val
            idx+=1
            last_line = line

        k=0
        dm_mo=[]
        bas_idx = 0

        # symmetric block first
        for num in self.nbpsy:
            arr = np.zeros([num-self.nfcpsy[bas_idx], num-self.nfcpsy[bas_idx]])
            if self.nfcpsy[bas_idx]> 0 :
                dm_mo=block_diag(dm_mo, occ*np.eye(self.nfcpsy[bas_idx]))
            for i in range(num-self.nfcpsy[bas_idx]):
                for j in range (i+1) :
                    arr[i, j]=float(sym_store[k])
                    if i!=j:
                        arr[j, i]=arr[i, j]
                        if symm_indices[k] in asymm_store:
                            asymm_val = asymm_store[symm_indices[k]]
                            arr[i,j] -= asymm_val
                            arr[j,i] += asymm_val
                    k+=1
            bas_idx +=1
            dm_mo=block_diag(dm_mo, arr)
        dm1_mo = np.delete(dm_mo, (0), axis=0)

        return dm1_mo


    def sdm_ao(self, i):
        fnameIN ='cid1trfl.FROMdrt1.state{}TOdrt1.state{}.iwfmt'.format(i, i)
        return self.mo_coeff @ self.dm_from_iwfmt(fnameIN, state_dm=True) @ self.mo_coeff.T

    def tdm_ao(self, iFROM, iTO):
        fnameIN ='cid1trfl.FROMdrt1.state{}TOdrt1.state{}.iwfmt'.format(iFROM, iTO)
        return self.mo_coeff @ self.dm_from_iwfmt(fnameIN, state_dm=False) @ self.mo_coeff.T

    def sdm_redc_ao(self,i):

        pass


    def H0(self, en_type):
        H0_arr=[]
        lines = [line for line in open('ciudgsm')]
        str_in = ['eci', 'eci+dv1', 'eci+dv2', 'eci+dv3', 'eci+pople']
        str_search = ['eci       =', 'eci+dv1   =', 'eci+dv2   =', 'eci+dv3   =', 'eci+pople =']

        if (en_type in str_in) == False:
            print('The following input strings are available:' ' \n {}'.format(str_in))
            exit
        else:
            str_lookup=str_search[str_in.index(en_type)]
            for line in lines:
                if 'convergence not reached' in line:
                    print ('\n''Convergence of all roots NOT reached. Throwing a WARNING here!!''\n')
                    continue
                elif str_lookup in line:
                    H0_arr = block_diag(H0_arr, float(line.split()[2]))
            self.H0_diag=np.delete(H0_arr, (0), axis=0)

        return self.H0_diag


if __name__ == '__main__':

    nstates=int(sys.argv[1])
    
    parser = colparser('molden_mo_mc.sp', 'tranls', nstates)
    nbas=parser.nbft
    print(parser.irrep_info())

    molden_dict = {"basis_file":"molden_mo_mc.sp",
    "molecule": "molden"}
    s = pyopencap.System(molden_dict)
    cap_dict = {"cap_type": "voronoi","r_cut": "3.00"}


    pc = pyopencap.CAP(s,cap_dict,nstates)
    pc.compute_ao_cap(cap_dict)
    W_ao = pc.get_ao_cap()
    ovlp = s.get_overlap_mat()
    H0 = parser.H0('eci')
    W = np.zeros((nstates,nstates))


    for i in range(0,nstates):
        for j in range(i,nstates):
            if i==j:
                print("{},{}".format(i+1,j+1))
                dm1_ao = parser.sdm_ao(i+1)
                print(np.trace(ovlp@dm1_ao))
            else:
                dm1_ao = parser.tdm_ao(i+1, j+1)

            W[i][j] = -1.0*np.trace(W_ao@dm1_ao)
            W[j][i] = -1.0*np.trace(W_ao@dm1_ao.T)
    print(pd.DataFrame(H0).to_string())
    print(pd.DataFrame(W).to_string())



    import matplotlib.pyplot as plt
    CAPH = CAPHamiltonian(H0=H0,W=W)
    ref_energy = np.min(H0)
    eta_list = np.linspace(0,1000,1001) * 1E-5
    CAPH.run_trajectory(eta_list)
    plt.plot(np.real(CAPH.energies_ev(ref_energy=ref_energy)),
        np.imag(CAPH.energies_ev(ref_energy=ref_energy)),'ro',label='Uncorrected Trajectory')
    plt.show()
