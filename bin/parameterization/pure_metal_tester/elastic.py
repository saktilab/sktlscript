#!/usr/bin/env python3

import sys
import pymatgen
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import scipy.constants as constants
import numpy as np
import os
import utils
import math
import copy

import time

class ElaticEvaluator():

    Ls_Dic={                       \
    '01':[ 1., 1., 1., 0., 0., 0.],\
    '02':[ 1., 0., 0., 0., 0., 0.],\
    '03':[ 0., 1., 0., 0., 0., 0.],\
    '04':[ 0., 0., 1., 0., 0., 0.],\
    '05':[ 0., 0., 0., 2., 0., 0.],\
    '06':[ 0., 0., 0., 0., 2., 0.],\
    '07':[ 0., 0., 0., 0., 0., 2.],\
    '08':[ 1., 1., 0., 0., 0., 0.],\
    '09':[ 1., 0., 1., 0., 0., 0.],\
    '10':[ 1., 0., 0., 2., 0., 0.],\
    '11':[ 1., 0., 0., 0., 2., 0.],\
    '12':[ 1., 0., 0., 0., 0., 2.],\
    '13':[ 0., 1., 1., 0., 0., 0.],\
    '14':[ 0., 1., 0., 2., 0., 0.],\
    '15':[ 0., 1., 0., 0., 2., 0.],\
    '16':[ 0., 1., 0., 0., 0., 2.],\
    '17':[ 0., 0., 1., 2., 0., 0.],\
    '18':[ 0., 0., 1., 0., 2., 0.],\
    '19':[ 0., 0., 1., 0., 0., 2.],\
    '20':[ 0., 0., 0., 2., 2., 0.],\
    '21':[ 0., 0., 0., 2., 0., 2.],\
    '22':[ 0., 0., 0., 0., 2., 2.],\
    '23':[ 0., 0., 0., 2., 2., 2.],\
    '24':[-1., .5, .5, 0., 0., 0.],\
    '25':[ .5,-1., .5, 0., 0., 0.],\
    '26':[ .5, .5,-1., 0., 0., 0.],\
    '27':[ 1.,-1., 0., 0., 0., 0.],\
    '28':[ 1.,-1., 0., 0., 0., 2.],\
    '29':[ 0., 1.,-1., 0., 0., 2.],\
    '30':[ .5, .5,-1., 0., 0., 2.],\
    '31':[ 1., 0., 0., 2., 2., 0.],\
    '32':[ 1., 1.,-1., 0., 0., 0.],\
    '33':[ 1., 1., 1.,-2.,-2.,-2.],\
    '34':[ .5, .5,-1., 2., 2., 2.],\
    '35':[ 0., 0., 0., 2., 2., 4.],\
    '36':[ 1., 2., 3., 4., 5., 6.],\
    '37':[-2., 1., 4.,-3., 6.,-5.],\
    '38':[ 3.,-5.,-1., 6., 2.,-4.],\
    '39':[-4.,-6., 5., 1.,-3., 2.],\
    '40':[ 5., 4., 6.,-2.,-1.,-3.],\
    '41':[-6., 3.,-2., 5.,-4., 1.]}


    Ls_str={                                     \
    '01':'(  eta,  eta,  eta,  0.0,  0.0,  0.0)',\
    '02':'(  eta,  0.0,  0.0,  0.0,  0.0,  0.0)',\
    '03':'(  0.0,  eta,  0.0,  0.0,  0.0,  0.0)',\
    '04':'(  0.0,  0.0,  eta,  0.0,  0.0,  0.0)',\
    '05':'(  0.0,  0.0,  0.0, 2eta,  0.0,  0.0)',\
    '06':'(  0.0,  0.0,  0.0,  0.0, 2eta,  0.0)',\
    '07':'(  0.0,  0.0,  0.0,  0.0,  0.0, 2eta)',\
    '08':'(  eta,  eta,  0.0,  0.0,  0.0,  0.0)',\
    '09':'(  eta,  0.0,  eta,  0.0,  0.0,  0.0)',\
    '10':'(  eta,  0.0,  0.0, 2eta,  0.0,  0.0)',\
    '11':'(  eta,  0.0,  0.0,  0.0, 2eta,  0.0)',\
    '12':'(  eta,  0.0,  0.0,  0.0,  0.0, 2eta)',\
    '13':'(  0.0,  eta,  eta,  0.0,  0.0,  0.0)',\
    '14':'(  0.0,  eta,  0.0, 2eta,  0.0,  0.0)',\
    '15':'(  0.0,  eta,  0.0,  0.0, 2eta,  0.0)',\
    '16':'(  0.0,  eta,  0.0,  0.0,  0.0, 2eta)',\
    '17':'(  0.0,  0.0,  eta, 2eta,  0.0,  0.0)',\
    '18':'(  0.0,  0.0,  eta,  0.0, 2eta,  0.0)',\
    '19':'(  0.0,  0.0,  eta,  0.0,  0.0, 2eta)',\
    '20':'(  0.0,  0.0,  0.0, 2eta, 2eta,  0.0)',\
    '21':'(  0.0,  0.0,  0.0, 2eta,  0.0, 2eta)',\
    '22':'(  0.0,  0.0,  0.0,  0.0, 2eta, 2eta)',\
    '23':'(  0.0,  0.0,  0.0, 2eta, 2eta, 2eta)',\
    '24':'( -eta,.5eta,.5eta,  0.0,  0.0,  0.0)',\
    '25':'(.5eta, -eta,.5eta,  0.0,  0.0,  0.0)',\
    '26':'(.5eta,.5eta, -eta,  0.0,  0.0,  0.0)',\
    '27':'(  eta, -eta,  0.0,  0.0,  0.0,  0.0)',\
    '28':'(  eta, -eta,  0.0,  0.0,  0.0, 2eta)',\
    '29':'(  0.0,  eta, -eta,  0.0,  0.0, 2eta)',\
    '30':'(.5eta,.5eta, -eta,  0.0,  0.0, 2eta)',\
    '31':'(  eta,  0.0,  0.0, 2eta, 2eta,  0.0)',\
    '32':'(  eta,  eta, -eta,  0.0,  0.0,  0.0)',\
    '33':'(  eta,  eta,  eta,-2eta,-2eta,-2eta)',\
    '34':'(.5eta,.5eta, -eta, 2eta, 2eta, 2eta)',\
    '35':'(  0.0,  0.0,  0.0, 2eta, 2eta, 4eta)',\
    '36':'( 1eta, 2eta, 3eta, 4eta, 5eta, 6eta)',\
    '37':'(-2eta, 1eta, 4eta,-3eta, 6eta,-5eta)',\
    '38':'( 3eta,-5eta,-1eta, 6eta, 2eta,-4eta)',\
    '39':'(-4eta,-6eta, 5eta, 1eta,-3eta, 2eta)',\
    '40':'( 5eta, 4eta, 6eta,-2eta,-1eta,-3eta)',\
    '41':'(-6eta, 3eta,-2eta, 5eta,-4eta, 1eta)'}

    LC_Dic = {              \
    'CI' :'Cubic I'        ,\
    'CII':'Cubic II'       ,\
    'HI' :'Hexagonal I'    ,\
    'HII':'Hexagonal II'   ,\
    'RI' :'Rhombohedral I' ,\
    'RII':'Rhombohedral II',\
    'TI' :'Tetragonal I'   ,\
    'TII':'Tetragonal II'  ,\
    'O'  :'Orthorhombic'   ,\
    'M'  :'Monoclinic'     ,\
    'N'  :'Triclinic'} 

    head = {                                                                     \
'CI':'\
    for, space group-number between 207 and 230, Cubic I structure.        \n\n\
               C11     C12     C12      0       0       0                  \n\
               C12     C11     C12      0       0       0                  \n\
               C12     C12     C11      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0      C44                 \n',\
'CII':'\
    for, space group-number between 195 and 206, Cubic II structure.       \n\n\
               C11     C12     C12      0       0       0                  \n\
               C12     C11     C12      0       0       0                  \n\
               C12     C12     C11      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0      C44                 \n',\
'HI':'\
    for, space group-number between 177 and 194, Hexagonal I structure.    \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C11     C13      0       0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0   (C11-C12)/2            \n',\
'HII':'\
    for, space group-number between 168 and 176, Hexagonal II structure.   \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C11     C13      0       0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0   (C11-C12)/2            \n',\
'RI':'\
    for, space group-number between 149 and 167, Rhombohedral I structure. \n\n\
               C11     C12     C13     C14      0       0                  \n\
               C12     C11     C13    -C14      0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
               C14    -C14      0      C44      0       0                  \n\
                0       0       0       0      C44     C14                 \n\
                0       0       0       0      C14  (C11-C12)/2            \n',\
'RII':'\
    for, space group-number between 143 and 148, Rhombohedral II structure.\n\n\
               C11     C12     C13     C14     C15      0                  \n\
               C12     C11     C13    -C14    -C15      0                  \n\
               C13     C13     C33      0       0       0                  \n\
               C14    -C14      0      C44      0     -C15                 \n\
               C15    -C15      0       0      C44     C14                 \n\
                0       0       0     -C15     C14  (C11-C12)/2            \n',\
'TI':'\
    for, space group-number between 89 and 142, Tetragonal I structure.    \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C11     C13      0       0       0                  \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
                0       0       0       0       0      C66                 \n',\
'TII':'\
    for, space group-number between 75 and 88, Tetragonal II structure.    \n\n\
               C11     C12     C13      0       0      C16                 \n\
               C12     C11     C13      0       0     -C16                 \n\
               C13     C13     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C44      0                  \n\
               C16    -C16      0       0       0      C66                 \n',\
'O':'\
    for, space group-number between 16 and 74, Orthorhombic structure.     \n\n\
               C11     C12     C13      0       0       0                  \n\
               C12     C22     C23      0       0       0                  \n\
               C13     C23     C33      0       0       0                  \n\
                0       0       0      C44      0       0                  \n\
                0       0       0       0      C55      0                  \n\
                0       0       0       0       0      C66                 \n',\
'M':'\
    for, space group-number between 3 and 15, Monoclinic structure.        \n\n\
               C11     C12     C13      0       0      C16                 \n\
               C12     C22     C23      0       0      C26                 \n\
               C13     C23     C33      0       0      C36                 \n\
                0       0       0      C44     C45      0                  \n\
                0       0       0      C45     C55      0                  \n\
               C16     C26     C36      0       0      C66                 \n',\
'N':'\
    for, space group-number between 1 and 2, Triclinic structure.          \n\n\
               C11     C12     C13     C14      C15    C16                 \n\
               C12     C22     C23     C24      C25    C26                 \n\
               C13     C23     C33     C34      C35    C36                 \n\
               C14     C24     C34     C44      C45    C46                 \n\
               C15     C25     C35     C45      C55    C56                 \n\
               C16     C26     C36     C46      C56    C66                 \n'}

    ordr = 2
    mthd = 'Energy'
    mdr = 0.05
    NoP = 21
    V0 = None
    LC = None
    ECs = None

    def make_elastic_structures(self, struct):
        sgroup = SpacegroupAnalyzer(struct)
        SGN = sgroup.get_space_group_number()
        SGN_explanation=sgroup.get_space_group_symbol()   

        if (1 <= SGN and SGN <= 2):      # Triclinic
            self.LC = 'N'
            self.ECs = 21
        elif(3 <= SGN and SGN <= 15):    # Monoclinic
            self.LC = 'M'
            self.ECs = 13
        elif(16 <= SGN and SGN <= 74):   # Orthorhombic
            self.LC = 'O'
            self.ECs =  9
        elif(75 <= SGN and SGN <= 88):   # Tetragonal II
            self.LC = 'TII'
            self.ECs =  7
        elif(89 <= SGN and SGN <= 142):  # Tetragonal I
            self.LC = 'TI'
            self.ECs =  6       
        elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
            self.LC = 'RII'
            self.ECs =  7
        elif(149 <= SGN and SGN <= 167): # Rhombohedral I
            self.LC = 'RI'
            self.ECs =  6
        elif(168 <= SGN and SGN <= 176): # Hexagonal II
            self.LC = 'HII'
            self.ECs =  5
        elif(177 <= SGN and SGN <= 194): # Hexagonal I
            self.LC = 'HI'
            self.ECs =  5
        elif(195 <= SGN and SGN <= 206): # Cubic II
            self.LC = 'CII'
            self.ECs =  3
        elif(207 <= SGN and SGN <= 230): # Cubic I
            self.LC = 'CI'
            self.ECs =  3
        else: sys.exit('\n.... Oops ERROR: WRONG Space-Group Number !?!?!?    \n')

        order = 'second'
        print ('\n     '+ SGN_explanation +'\
            \n     '+ self.LC_Dic[self.LC] +' structure in the Laue classification.\
            \n     This structure has '+ str(self.ECs) +' independent '+ order +'-order elastic constants.')
        #--------------------------------------------------------------------------------------------------

        
        mdr = round(self.mdr, 3)
        print ('     The maximum Lagrangian strain is '+ str(mdr) + '\n')
        #--------------------------------------------------------------------------------------------------

        
        
        print ('     The number of the distorted structures is '+ str(self.NoP) + '\n')

        ptn = int((self.NoP-1)/2)

        interval = 0.0001

        if (mdr/ptn <= interval):
            sys.exit('.... Oops ERROR: The interval of the strain values is < '+ str(interval) +\
                '\n                 Choose a larger maximum Lagrangian strain'\
                '\n                 or a less number of distorted structures.\n')
        #--------------------------------------------------------------------------------------------------

        #%%-- Reading the scale, stretch, and base vectors from the input file and calculate the volume -%%
        scale = 1.0
        stretch=[1.,1.,1.]
        a0 = constants.physical_constants['Bohr radius'][0]/constants.angstrom
        
        bv = struct.lattice.matrix/a0

        M_old= np.array(bv)
        D    = np.linalg.det(M_old)
        self.V0   = abs(stretch[0]*stretch[1]*stretch[2]*scale**3*D)
        #--------------------------------------------------------------------------------------------------
        
        if (self.LC == 'CI' or \
            self.LC == 'CII'):
            Lag_strain_list = ['01','08','23']
        if (self.LC == 'HI' or \
            self.LC == 'HII'):
            Lag_strain_list = ['01','26','04','03','17']
        if (self.LC == 'RI'):
            Lag_strain_list = ['01','08','04','02','05','10']
        if (self.LC == 'RII'):
            Lag_strain_list = ['01','08','04','02','05','10','11']
        if (self.LC == 'TI'):
            Lag_strain_list = ['01','26','27','04','05','07']
        if (self.LC == 'TII'):
            Lag_strain_list = ['01','26','27','28','04','05','07']
        if (self.LC == 'O'):
            Lag_strain_list = ['01','26','25','27','03','04','05','06','07']
        if (self.LC == 'M'):
            Lag_strain_list = ['01','25','24','28','29','27','20','12','03','04','05','06','07']
        if (self.LC == 'N'):
            Lag_strain_list = ['02','03','04','05','06','07','08','09','10','11',\
                            '12','13','14','15','16','17','18','19','20','21','22']

        
        cont1= 0

        distorted_structures = {}


        for i in Lag_strain_list:
            Ls_list= self.Ls_Dic[i]

            cont1  = cont1 + 1
            if (cont1 < 10):
                Dstn = 'Dst0'+str(cont1)
            else:
                Dstn = 'Dst' +str(cont1)
            
            distorted_structures[Dstn] = []
            
            cont2 = 0
            for s in range(-ptn, ptn+1):
                r = mdr*s/ptn
                if (s==0):
                    r = 0.0001
                    

                Ls = np.zeros(6)
                for i in range(6):
                    Ls[i] = Ls_list[i]
                Lv = r*Ls
                # Lagrangian strain to physical strain (eta = eps + 0.5*eps*esp) --------------------------
                eta_matrix      = np.zeros((3,3))

                eta_matrix[0,0] = Lv[0]
                eta_matrix[0,1] = Lv[5]/2.
                eta_matrix[0,2] = Lv[4]/2.
                
                eta_matrix[1,0] = Lv[5]/2.
                eta_matrix[1,1] = Lv[1]
                eta_matrix[1,2] = Lv[3]/2.

                eta_matrix[2,0] = Lv[4]/2.
                eta_matrix[2,1] = Lv[3]/2.
                eta_matrix[2,2] = Lv[2]

                norm       = 1.0

                eps_matrix = eta_matrix
                if (np.linalg.norm(eta_matrix) > 0.7):
                    sys.exit('\n.... Oops ERROR: Too large deformation!\n') 

                while( norm > 1.e-10 ):
                    x          = eta_matrix - np.dot(eps_matrix, eps_matrix)/2.
                    norm       = np.linalg.norm(x - eps_matrix)      
                    eps_matrix = x

        #--------------------------------------------------------------------------------------------------
                i_matrix   = np.array([[1., 0., 0.],
                                    [0., 1., 0.], 
                                    [0., 0., 1.]])
                def_matrix = i_matrix + eps_matrix
                M_new      = np.dot(M_old, def_matrix)*a0
        
                new_latt = pymatgen.Lattice(M_new) 
                # det    = np.linalg.det(M_new)
                


                new_str = struct.copy()
                new_str.modify_lattice(new_latt)

                distorted_structures[Dstn].append( (r, new_str) )
        
        return distorted_structures

    def sortlist(self, lst1, lst2):
        temp = copy.copy(lst1)

        lst3 = []
        lst4 = []

        temp.sort()

        for i in range(len(lst1)):
            lst3.append(lst1[lst1.index(temp[i])])
            lst4.append(lst2[lst1.index(temp[i])])

        return lst3, lst4

    def analyze_elastic_constants(self, evaluated_res):
        _e     = 1.602176565e-19              # elementary charge
        Bohr   = 5.291772086e-11              # a.u. to meter
        Ryd2eV = 13.605698066                 # Ryd to eV
        Hartree2eV = 27.21138602              # Hartree to eV
        ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
        CONV = ToGPa * math.factorial(self.ordr)*2.

        for i in range(1, self.ECs+1):
            if (i<10):
                Dstn = 'Dst0'+str(i)
            else:
                Dstn = 'Dst' +str(i)
                
                if (Dstn not in evaluated_res):
                    raise ValueError('Wrong')

            
            for j in [6, 4]: #range(self.ordr+4, self.ordr-1, -2):
                strain =  [x[0] for x in evaluated_res[Dstn]]
                energy =  [x[1] for x in evaluated_res[Dstn]]
                strain, energy = self.sortlist(strain, energy)
                strain0 = copy.copy(strain)
                energy0 = copy.copy(energy)
                # print('Deformation: {}, Fitting order: {}'.format(Dstn, j))        
                # print('# Max. eta    SUM(Cij)')
                while (len(strain) > j): 
                    emax  = max(strain)
                    emin  = min(strain)
                    emax  = max(abs(emin),abs(emax))
                    coeffs= np.polyfit(strain, energy, j)
                    
                    Cij  = coeffs[j-2]*CONV/self.V0         # in GPa unit 
                    

                    # print('%13.10f'%emax, '%18.6f'%Cij)

                    if (abs(strain[0]+emax) < 1.e-7):
                        strain.pop(0); energy.pop(0)
                    if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                        strain.pop()
                        energy.pop()

                #--- Cross-Validation error calculations --------------------------------------------------
                strain = copy.copy(strain0)
                energy = copy.copy(energy0)

                
                # print('# Max. eta    Cross-Validation Error   ')
                while (len(strain) > j+1): 
                    emax = max(strain)
                    emin = min(strain)
                    emax = max(abs(emin),abs(emax))

                    S = 0
                    for k in range(len(strain)):
                        Y      = energy[k]
                        etatmp = []
                        enetmp = []

                        for l in range(len(strain)):
                            if (l==k): pass
                            else:            
                                etatmp.append(strain[l])
                                enetmp.append(energy[l])

                        Yfit = np.polyval(np.polyfit(etatmp,enetmp, j), strain[k])
                        S    = S + (Yfit-Y)**2

                    CV = math.sqrt(S/len(strain))
                    # print('%13.10f'%emax, CV)

                    if (abs(strain[0]+emax) < 1.e-7):
                        strain.pop(0)
                        energy.pop(0)
                    if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                        strain.pop()
                        energy.pop()
        A2 = []
        for i in range(1, self.ECs+1):
            if (i<10):
                Dstn = 'Dst0'+str(i)
            else:
                Dstn = 'Dst' +str(i)

            mdri   = 0.03
            ordri  = 6

            strain = []
            energy = []

            for stri, energy_ in evaluated_res[Dstn]:
                if (abs(stri) <= mdri):
                    strain.append(stri)
                    energy.append(energy_)
            # strain =  [x[0] for x in evaluated_res[Dstn]]
            # energy =  [x[1] for x in evaluated_res[Dstn]]


            coeffs = np.polyfit(strain, energy, ordri)
            A2.append(coeffs[ordri-2])

        A2 = np.array(A2)
        if (len(A2) != self.ECs):
            sys.exit('\n.... Oops ERROR: The number of data in the "ElaStic_2nd.in" is NOT equal to ' + \
            str(self.ECs)+'\n')

        C = np.zeros((6,6))

        #%!%!%--- Cubic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (self.LC == 'CI' or \
            self.LC == 'CII'):
            C[0,0] =-2.*(A2[0]-3.*A2[1])/3.
            C[1,1] = C[0,0]
            C[2,2] = C[0,0]
            C[3,3] = A2[2]/6.
            C[4,4] = C[3,3]
            C[5,5] = C[3,3]
            C[0,1] = (2.*A2[0]-3.*A2[1])/3.
            C[0,2] = C[0,1]
            C[1,2] = C[0,1]
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Hexagonal structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (self.LC == 'HI' or \
            self.LC == 'HII'):
            C[0,0] = 2.*A2[3]
            C[0,1] = 2./3.*A2[0] + 4./3.*A2[1] - 2.*A2[2] - 2.*A2[3]
            C[0,2] = 1./6.*A2[0] - 2./3.*A2[1] + 0.5*A2[2]
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[2,2] = 2.*A2[2]
            C[3,3] =-0.5*A2[2] + 0.5*A2[4]
            C[4,4] = C[3,3]
            C[5,5] = .5*(C[0,0] - C[0,1])
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Rhombohedral I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        if (self.LC == 'RI'):
            C[0,0] = 2.*A2[3]
            C[0,1] = A2[1]- 2.*A2[3]
            C[0,2] = .5*( A2[0] - A2[1] - A2[2])
            C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[1,3] =-C[0,3]
            C[2,2] = 2.*A2[2]
            C[3,3] = .5*A2[4]
            C[4,4] = C[3,3]
            C[4,5] = C[0,3]
            C[5,5] = .5*(C[0,0] - C[0,1])
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Rhombohedral II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (self.LC == 'RII'):
            C[0,0] = 2.*A2[3]
            C[0,1] = A2[1]- 2.*A2[3]
            C[0,2] = .5*( A2[0] - A2[1] - A2[2])
            C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
            C[0,4] = .5*(-A2[3] - A2[4] + A2[6])
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[1,3] =-C[0,3]
            C[1,4] =-C[0,4]    
            C[2,2] = 2.*A2[2]
            C[3,3] = .5*A2[4]
            C[3,5] =-C[0,4]
            C[4,4] = C[3,3]
            C[4,5] = C[0,3]
            C[5,5] = .5*(C[0,0] - C[0,1])
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Tetragonal I structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        if (self.LC == 'TI'):
            C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[3]
            C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[3]
            C[0,2] = A2[0]/6.-2.*A2[1]/3.+.5*A2[3]
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[2,2] = 2.*A2[3]
            C[3,3] = .5*A2[4]
            C[4,4] = C[3,3]
            C[5,5] = .5*A2[5]
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Tetragonal II structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (self.LC == 'TII'):
            C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[4]
            C[1,1] = C[0,0]
            C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[4]
            C[0,2] = A2[0]/6.-(2./3.)*A2[1]+.5*A2[4]
            C[0,5] = (-A2[2]+A2[3]-A2[6])/4.
            C[1,2] = C[0,2]
            C[1,5] =-C[0,5]
            C[2,2] = 2.*A2[4]
            C[3,3] = .5*A2[5]
            C[4,4] = C[3,3]
            C[5,5] = .5*A2[6]
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Orthorhombic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        if (self.LC == 'O'):
            C[0,0] = 2.*A2[0]/3.+4.*A2[1]/3.+A2[3]-2.*A2[4]-2.*A2[5]
            C[0,1] = 1.*A2[0]/3.+2.*A2[1]/3.-.5*A2[3]-A2[5]
            C[0,2] = 1.*A2[0]/3.-2.*A2[1]/3.+4.*A2[2]/3.-.5*A2[3]-A2[4]
            C[1,1] = 2.*A2[4]
            C[1,2] =-2.*A2[1]/3.-4.*A2[2]/3.+.5*A2[3]+A2[4]+A2[5]
            C[2,2] = 2.*A2[5]
            C[3,3] = .5*A2[6]
            C[4,4] = .5*A2[7]
            C[5,5] = .5*A2[8]
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Monoclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        if (self.LC == 'M'):
            C[0,0] = 2.*A2[0]/3.+8.*(A2[1]+A2[2])/3.-2.*(A2[5]+A2[8]+A2[9])
            C[0,1] = A2[0]/3.+4.*(A2[1]+A2[2])/3.-2.*A2[5]-A2[9]
            C[0,2] =(A2[0]-4.*A2[2])/3.+A2[5]-A2[8]
            C[0,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.+.5*(A2[5]+A2[7]+A2[8]+A2[9]-A2[12])
            C[1,1] = 2.*A2[8]
            C[1,2] =-4.*(2.*A2[1]+A2[2])/3.+2.*A2[5]+A2[8]+A2[9]
            C[1,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.-.5*A2[3]+A2[5]+.5*(A2[7]+A2[8]+A2[9])
            C[2,2] = 2.*A2[9]
            C[2,5] =-1.*A2[0]/6.+2.*A2[1]/3.-.5*(A2[3]+A2[4]-A2[7]-A2[8]-A2[9]-A2[12])
            C[3,3] = .5*A2[10]
            C[3,4] = .25*(A2[6]-A2[10]-A2[11])
            C[4,4] = .5*A2[11]
            C[5,5] = .5*A2[12]
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Triclinic structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (self.LC == 'N'):
            C[0,0] = 2.*A2[0]
            C[0,1] = 1.*(-A2[0]-A2[1]+A2[6])
            C[0,2] = 1.*(-A2[0]-A2[2]+A2[7])
            C[0,3] = .5*(-A2[0]-A2[3]+A2[8]) 
            C[0,4] = .5*(-A2[0]+A2[9]-A2[4])
            C[0,5] = .5*(-A2[0]+A2[10]-A2[5])
            C[1,1] = 2.*A2[1]
            C[1,2] = 1.*(A2[11]-A2[1]-A2[2])
            C[1,3] = .5*(A2[12]-A2[1]-A2[3])
            C[1,4] = .5*(A2[13]-A2[1]-A2[4])
            C[1,5] = .5*(A2[14]-A2[1]-A2[5])
            C[2,2] = 2.*A2[2] 
            C[2,3] = .5*(A2[15]-A2[2]-A2[3])
            C[2,4] = .5*(A2[16]-A2[2]-A2[4])
            C[2,5] = .5*(A2[17]-A2[2]-A2[5])
            C[3,3] = .5*A2[3]
            C[3,4] = .25*(A2[18]-A2[3]-A2[4])
            C[3,5] = .25*(A2[19]-A2[3]-A2[5])
            C[4,4] = .5*A2[4]
            C[4,5] = .25*(A2[20]-A2[4]-A2[5])
            C[5,5] = .5*A2[5]
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Calculating the elastic moduli ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        CONV = ToGPa * 2.    
        

        for i in range(5):
            for j in range(i+1,6):
                C[j,i] = C[i,j] 

        C = C * CONV/self.V0

        BV = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9
        GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
        EV = (9*BV*GV)/(3*BV+GV)
        nuV= (1.5*BV-GV)/(3*BV+GV)
        S  = np.linalg.inv(C)
        BR = 1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]))
        GR =15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
        ER = (9*BR*GR)/(3*BR+GR)
        nuR= (1.5*BR-GR)/(3*BR+GR)
        BH = 0.50*(BV+BR)
        GH = 0.50*(GV+GR)
        EH = (9.*BH*GH)/(3.*BH+GH)
        nuH= (1.5*BH-GH)/(3.*BH+GH)
        AVR= 100.*(GV-GR)/(GV+GR)
        #--------------------------------------------------------------------------------------------------
        lineseparator = '%'*80

        #%!%!%--- Writing the output file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        print('    Today is '+ time.asctime() +                                           '\n'\
                                                                                                '\n'\
                    '    Symmetry of the second-order elastic constant matrix in Voigt notation. \n'\
                        + self.head[self.LC] +                                                             '\n'\
                    '    Elastic constant (stiffness) matrix in GPa:                             \n')

        for i in range(0,6):
            print('', end=' ')
            for j in range(0,6):
                print('%11.1f'%(C[i,j]), end=' ')
            print()

        # print('\n\n    Elastic compliance matrix in 1/GPa: \n')

        # for i in range(0,6):
        #     print('', end=' ')
        #     for j in range(0,6):
        #         print('%11.5f'%(S[i,j]), end=' ')
        #     print()

        print('\n'+ lineseparator +'\n')

        print('    Voigt bulk  modulus, B_V = {0}  GPa'.format('%8.2f'%(BV)))
        print('    Voigt shear modulus, G_V = {0}  GPa'.format('%8.2f'%(GV)) + '\n')

        print('    Reuss bulk  modulus, B_R = {0}  GPa'.format('%8.2f'%(BR)))
        print('    Reuss shear modulus, G_R = {0}  GPa'.format('%8.2f'%(GR)) + '\n')

        print('    Hill bulk  modulus,  B_H = {0}  GPa'.format('%8.2f'%(BH)))
        print('    Hill shear modulus,  G_H = {0}  GPa'.format('%8.2f'%(GH)))

        print('\n'+ lineseparator +'\n')

        print('    Voigt Young modulus,  E_V = {0}  GPa'.format('%8.2f'%(EV)))
        print('    Voigt Poisson ratio, nu_V = {0}'     .format('%8.2f'%(nuV)) + '\n')

        print('    Reuss Young modulus,  E_R = {0}  GPa'.format('%8.2f'%(ER)))
        print('    Reuss Poisson ratio, nu_R = {0}'     .format('%8.2f'%(nuR)) + '\n')

        print('    Hill Young modulus,   E_H = {0}  GPa'.format('%8.2f'%(EH)))
        print('    Hill Poisson ratio,  nu_H = {0}'     .format('%8.2f'%(nuH)))

        print('\n'+ lineseparator +'\n')

        print('    Elastic Anisotropy in polycrystalline, AVR = {0} %'.format('%8.3f'%(AVR)))

        print('\n'+ lineseparator +'\n')

        print('    Eigenvalues of elastic constant (stiffness) matrix:   \n')

        eigval=np.linalg.eig(C)
        for i in range(6):
            print('%16.1f' % float(eigval[0][i]))

        print('\n'+ lineseparator +'\n')

        res = {}
        res['B_V'] = BV
        res['G_V'] = GV
        res['B_R'] = BR
        res['G_R'] = GR
        res['B_H'] = BH
        res['G_H'] = GH
        res['E_V'] = EV
        res['nu_V'] = nuV
        res['E_R'] = ER
        res['nu_R'] = nuR
        res['E_H'] = EH
        res['nu_H'] = nuH
        res['AVR'] = AVR

        return res    

    def compute_elastic_constants(self, root_path, struct, skfile_path, temp_input_path, dftbplus_path):
        running_path = os.path.join(root_path, 'elastic')
        os.makedirs(running_path, exist_ok=True)
        all_distorted_structures = self.make_elastic_structures(struct)
        res = {}
        for dist, strain_structs in all_distorted_structures.items():
            print('Computing deformation {}'.format(dist))
            fout = open(os.path.join(running_path, '{}_Energy.dat'.format(dist)), 'w')
            res[dist] = []
            for ind, strain_struct in enumerate(strain_structs):
                strain, struct = strain_struct
                print('   Evaluate {}/{}...'.format(ind+1, len(strain_structs)))
                temp_path = os.path.join(running_path, '{}_{:02d}'.format(dist, ind))
                l_res = utils.evaluate_dftb(struct, temp_path, skfile_path=skfile_path, temp_input_path=temp_input_path, dftbplus_path=dftbplus_path, geom_opt=False)
                # print(ind, strain, l_res['energy_final'])
                print('{:20.12f} {:20.12f}'.format(strain, l_res['energy_final']), file=fout)
                res[dist].append( (strain, l_res['energy_final']))
            fout.close()
        elastic_res = self.analyze_elastic_constants(res)
        return elastic_res