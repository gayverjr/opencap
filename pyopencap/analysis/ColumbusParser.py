import os, sys
import re
import numpy as np
import pandas as pd
pd.set_option("display.precision", 24)
from pathlib2 import Path
from scipy.linalg import block_diag
from pandas import DataFrame as df


class ColumbusDM:

    def __init__(self, nstates):
        
        self.nstates = nstates


        rfile = open('WORK/tranls', 'r')
        while True:
            try:
                line=next(rfile)
            except StopIteration:
                break
            if 'transformation information:' in line:
                line=next(rfile)
                line=next(rfile)
                self.nbft = re.findall(r'\d+', line)[0]
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

    def irrep_info(self):
        pformat = r' Total number of basis functions: {}' \
        ' \n Symmetry Label: {}' \
        ' \n NBPSY: {}'' \n NMPSY: {}' ' \n NFCPSY: {} \n'
        return pformat.format(self.nbft, self.slabel, self.nbpsy, self.nmpsy, self.nfcpsy)

    
    def dm_from_iwfmt(self, fnameIN, ftype):
        if ftype=='tdm':occ = 0.0
        elif ftype=='sdm':occ = 2.0
        Nw=[]
        dm_store=[]
        syms = []
        sym = 1.0

        lines = [line for line in open(fnameIN)]
        for i, line in enumerate(lines):
            if ' 0.000000000000E+00' in line and len(line.split())==1:
                line = lines[i-1]
                (numX,lab1,ibvtyp,itypea,itypeb,ifmt,last,nipv) = [int(word) for word in line.split()]
                if itypea == 0 & itypeb == 7:
                    sym = 1.0
                elif itypea==2 & itypeb ==9 :
                    sym = -1.0
            elif len(line.split())==3:
                if line.split()[0].startswith( 'tout' ):
                    continue
                else:
                    dm_store.append(line.split()[0])
                    syms.append(sym)


        bas_index=0
        k=0
        for num in self.nbpsy:
            arr = np.zeros([num-self.nfcpsy[bas_index], num-self.nfcpsy[bas_index]])
            if self.nfcpsy[bas_index]> 0 :
                Nw=block_diag(Nw, occ*np.eye(self.nfcpsy[bas_index]))
            for i in range(num-self.nfcpsy[bas_index]):
                for j in range (i+1) :
                    arr[i, j]=float(dm_store[k]) *syms[k]
                    arr[j, i]=arr[i, j]
                    k+=1
            bas_index +=1
            Nw=block_diag(Nw, arr)
        Nw = np.delete(Nw, (0), axis=0)

        return Nw

    def sdm_mat(self, i):
        fnameIN ='cid1trfl.FROMdrt1.state{}TOdrt1.state{}.iwfmt'.format(i, i)
        sdm_mo = self.dm_from_iwfmt(fnameIN, 'sdm')
        return sdm_mo

    def tdm_mat(self, iFROM, iTO):
        fnameIN ='cid1trfl.FROMdrt1.state{}TOdrt1.state{}.iwfmt'.format(iFROM, iTO)
        tdm_mo = self.dm_from_iwfmt(fnameIN, 'tdm')
        return tdm_mo

    def dm_filegen(self, fnameIN, fnameOUT, ftype ):
        dm_mo = self.dm_from_iwfmt(fnameIN, ftype)

        fname = open(fnameOUT, 'w')
        print(df(dm_mo).to_string(index=False,header=False), file=fname)
        fname.close()

        return 'Generated: %s'%fnameOUT


if __name__ == '__main__':

    nstates=int(sys.argv[1])

    density_file=ColumbusDM(nstates)
    nbas=density_file.nbft
    print(density_file.irrep_info())
    

 #   '''
    #Print MO file
    mofile=open('molden_mo_mc.sp')
    lines = [line for line in mofile]
    for line in lines:
        if line.startswith('[MO]'):
            prnt_index = (lines.index(line))+1
            break
    mo_in=[]
    for i in range(prnt_index, len(lines)):
        mo_in.append(lines[i].split()[1])
    mo_save = np.reshape(mo_in, (nbas, nbas+4)).T
    mo_new=mo_save[4:nbas+4, :]
    fname=open('mocoef-file.dat', 'w')
    
    print(df(mo_new).to_string(index=False, header=False), file=fname)
    fname.close()
#    '''

#    '''
    #Print DM files from CI calculation
    for isource in range (1, nstates+1):
        for itarget in range (isource+1, nstates+1):
            fnameIN = 'cid1trfl.FROMdrt1.state{}TOdrt1.state{}.iwfmt'.format(isource, itarget)
            if os.path.isfile(fnameIN)==True :   #if bra and ket states are same, i-->i tdm files won't be there!
                fnameOUT = 'tdm.s{}TOs{}.dat'.format(isource, itarget)
                print(density_file.dm_filegen(fnameIN, fnameOUT, 'tdm'))
        

    for i in range (1, nstates+1):
            fnameIN ='cid1trfl.FROMdrt1.state{}TOdrt1.state{}.iwfmt'.format(i, i)
            fnameOUT = 'sdm{}.dat'.format(i)
            print(density_file.dm_filegen(fnameIN, fnameOUT, 'sdm'))
#    '''