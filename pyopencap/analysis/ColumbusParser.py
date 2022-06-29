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

import re
import numpy as np
from scipy.linalg import block_diag
import os

class MO:
    '''
    A class that is internally needed to read file containing MO coefficients.
    
    Used in read_molden function.

    ...

    Attributes
    ----------
    line : str
        Line string read from molden file
    total_idx : int
        Total no of line strings.

    '''
    def __init__(self,line,total_idx):
        '''
        Initializes the MO class

        Parameters
        ----------
        line : str
            Line string
        total_idx : int
            Total no of such line strings.
        '''

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


class colparser():
    '''
    A class that parses COLUMBUS electronic structure package generated files
    to generate state density and transition density matrices for projected-CAP
    calculation on MR-CI level.
    '''

    def _set_mo_coeff(self,ordering=None):
        '''
        An internal function to generate mo coffiecients in proper ordering.

        Returns
        -------
        MO coefficients in proper ordering to that of pyopencap internal structure
        '''

        if ordering is not None:
            raise NotImplementedError("Custom ordering NYI.")
        else:
            mo_coeff = np.zeros((self.nbft,self.nbft))
            for i, mo in enumerate(self.mos):
                mo_coeff[:,i] = mo.coeffs
            self.mo_coeff = mo_coeff


    def _read_molden(self,molden_file):
        '''
        An internal function that parses standard molden file.

        Parameters
        ----------
        molden_file : str
            molden MO filename (generated in MOLDEN/ folder in COLUMBUS calculation Directory
            if 'molden' keyword is invoked in control.run runfile.))

        Returns
        -------
        mos : np.ndarray
            MO coefficients.
        '''

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


    def __init__(self, molden_file, tranls):
        '''
        Initializes the colparser class

        Parameters
        ----------
        molden_file: str
            molden MO filename (generated in MOLDEN/ folder in COLUMBUS calculation Directory)
        tranls: str
            Path to tranls file generated by Columbus in WORK directory

        '''
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

        self.mos = self._read_molden(molden_file)
        self._set_mo_coeff()


    def mo_summary(self):
        '''
        Prints information about active space and symmetries of MOs.

        Returns
        -------
        str
            Total number of basis functions (nbft)
            Number of basis functions in each symmetry block (NBPSY).
            Number of orbitals in each of the symmetry blocks.(NMPSY<= NBSPSY)
            Number of frozen orbitals in each of the symmetry blocks (NFCSY ).
            Character labels for the symmetry blocks (SLABEL).
        '''

        pformat = r' Total number of basis functions: {}' \
        ' \n Symmetry Label: {}' \
        ' \n NBPSY: {}'' \n NMPSY: {}' ' \n NFCPSY: {} \n'
        return pformat.format(self.nbft, self.slabel, self.nbpsy, self.nmpsy, self.nfcpsy)

    
    def _dm_from_iwfmt(self, finame, state_dm=False):
        '''
        An internal parser that reads CI densities from user genreated 'iwfmt' density files and
        formats them in proper order to generate density matrices.

        Parameters
        ----------
        finame : str 
            File name from where the density elements needs to be read.
        state_dm : bool
            To specify which kind of density file are being parsed 
            ('True' for State density, 'False' for Transition density).

        Returns 
        -------
        dm1_mo: np.ndarray
            1 particle density matrix in MO basis
        '''

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

    def sdm_ao(self, i, DRTn=1, data_dir='.', filename=None):
        '''
        Returns state density matrix in atomic orbital basis by parsing a Columbus cid1trfl.iwfmt file.

        Parameters
        ----------
        i : int
            State index
        DRTn : int, optional
            DRT index
        data_dir: str, optional
            Directory to search for .iwfmt file. Should not be used in conjunction with `filename` kwarg
        filename: str, optional
            Path to file to parse. If not specified, the filename is assumed to be cid1trfl.FROMdrt{drtFrom}.state{i}TOdrt{drtTo}.state{i}.iwfmt 
            in the current directory.

        Returns
        -------
        sdm: np.ndarray
            State density matrix in AO basis

        '''
        if filename is None:
            fnameIN =os.path.join(data_dir,'cid1trfl.FROMdrt{}.state{}TOdrt{}.state{}.iwfmt'.format(DRTn, i, DRTn, i))
        else:
            fnameIN = filename
        return self.mo_coeff @ self._dm_from_iwfmt(fnameIN, state_dm=True) @ self.mo_coeff.T

    def tdm_ao(self, iFROM, iTO, drtFrom=1, drtTo=1, data_dir='.', filename=None):
        '''
        Returns transition density matrix in atomic orbital basis by parsing a Columbus cid1trfl.iwfmt file.

        Parameters
        ----------
        iFROM, iTO : int
            Initial state index and final state indices respectively.
        drtFrom, drtTo : int, optional
            DRT indices
        data_dir: str, optional
            Directory to search for .iwfmt file. Should not be used in conjunction with `filename` kwarg
        filename: str, optional
            Path to file to parse. If not specified, the filename is assumed to be cid1trfl.FROMdrt{drtFrom}.state{iFrom}TOdrt{drtTo}.state{iTO}.iwfmt 
            in the current directory.

        Returns
        -------
        tdm: np.ndarray
            Transition density matrix in AO basis

        '''
        if filename is None:
            fnameIN =os.path.join(data_dir,'cid1trfl.FROMdrt{}.state{}TOdrt{}.state{}.iwfmt'.format(drtFrom, iFROM, drtTo, iTO))
        else:
            fnameIN = filename
        return self.mo_coeff @ self._dm_from_iwfmt(fnameIN, state_dm=False) @ self.mo_coeff.T

    def sdm_ao_cid1fl(self,i, DRTn):
        '''
        Read CI state densities from state density cid1fl*.iwfmt files.

        Currently NYI.

        Parameters
        ----------
        i : int
            State index
        DRTn : int
            DRT index

        Raises
        -------
        NotImplementedError

        Returns
        -------
        np.ndarray
            State density matrix in AO basis
        '''
        raise NotImplementedError()


    def _H0_parser(self, str_arg, filename):
        '''
        An internal parser that parses 'ciudgsm' file and is used by H0 function.
        
        Parameters
        ----------
        str_arg : str
            One of the following: {'eci', 'eci+dv1', 'eci+dv2', 'eci+dv3', 'eci+pople' }
        filename: str
            Path to file

        Returns
        --------
        H0_diag: np.ndarray
            Diagonal Hamiltonian with CI energies.
        '''

        H0_arr=[]
        lines = [line for line in open(filename)]
        str_in = ['eci', 'eci+dv1', 'eci+dv2', 'eci+dv3', 'eci+pople']
        str_search = ['eci       =', 'eci+dv1   =', 'eci+dv2   =', 'eci+dv3   =', 'eci+pople =']
        str_lookup=str_search[str_in.index(str_arg)]
        for line in lines:
            if str_lookup in line:
                H0_arr = block_diag(H0_arr, float(line.split()[2]))
        self.H0_diag=np.delete(H0_arr, (0), axis=0)
        return self.H0_diag


    def get_H0(self, correction_type='eci+pople', filename='ciudgsm'):
        '''
        Parses energies from a Columbus ciudgsm file. 

        Parameters
        ----------
        correction_type: str, optional
            One of  {'eci+pople', 'eci', 'eci+dv1', 'eci+dv2', 'eci+dv3' }. Default is 'eci+pople'
        filename: str, optional
            Path to Columbus ciudgsm file located in WORK directory. If unspecified, assumed to be './ciudgsm'.

        Notes
        ------
        See https://aip.scitation.org/doi/pdf/10.1063/1.5144267 for discussion of corrections.

        Returns
        --------
        H0_mat : np.ndarray
            Diagonal hamiltonian with CI energies.

        '''
        str_in = ['eci', 'eci+dv1', 'eci+dv2', 'eci+dv3', 'eci+pople']
        if ((correction_type in str_in) == False) or (correction_type==None):
            print("The following input args are available (Switching to default 'eci+pople'):"  "\n {}".format(str_in))
            self.H0_mat = self._H0_parser('eci+pople',filename)
        else:
            self.H0_mat = self._H0_parser(correction_type,filename)

        return self.H0_mat
