#!/usr/bin/env python3
import sys
from collections import Counter 
import itertools
import os
import tempfile
import shutil
import subprocess
import re


SKPATH_ROOT = '/users/solccp/skfiles_dc/'
DCDFTBK_PATH = '/users/ynishi/bin/dftbk.Jul00.x'
SKFILE_SUFFIX = 'skf'

parameter_set = sys.argv[1]

SKPATH = SKPATH_ROOT + '/' + parameter_set

if (not os.path.isdir(SKPATH)):
    print ('The first argument is not a valid path for the SK files')
    sys.exit(1)

template_input = sys.argv[2]

if (not os.path.exists(template_input)):
    print ('The second argument is not a valid template DC-DFTB-K input')
    sys.exit(1)

system = sys.argv[3]
input_file = sys.argv[4]
output_file = sys.argv[5]

Element_symbol = {
    1: "H",    2: "He",    3: "Li",    4: "Be",    5: "B",    6: "C",
    7: "N",    8: "O",    9: "F",    10: "Ne",    11: "Na",    12: "Mg",
    13: "Al",    14: "Si",    15: "P",    16: "S",    17: "Cl",    18: "Ar",
    19: "K",    20: "Ca",    21: "Sc",    22: "Ti",    23: "V",    24: "Cr",
    25: "Mn",    26: "Fe",    27: "Co",    28: "Ni",    29: "Cu",    30: "Zn",
    31: "Ga",    32: "Ge",    33: "As",    34: "Se",    35: "Br",    36: "Kr",
    37: "Rb",    38: "Sr",    39: "Y",    40: "Zr",    41: "Nb",    42: "Mo",
    43: "Tc",    44: "Ru",    45: "Rh",    46: "Pd",    47: "Ag",    48: "Cd",
    49: "In",    50: "Sn",    51: "Sb",    52: "Te",    53: "I",    54: "Xe",
    55: "Cs",    56: "Ba",    57: "La",    58: "Ce",    59: "Pr",    60: "Nd",
    61: "Pm",    62: "Sm",    63: "Eu",    64: "Gd",    65: "Tb",    66: "Dy",
    67: "Ho",    68: "Er",    69: "Tm",    70: "Yb",    71: "Lu",    72: "Hf",
    73: "Ta",    74: "W",    75: "Re",    76: "Os",    77: "Ir",    78: "Pt",
    79: "Au",    80: "Hg",    81: "Tl",    82: "Pb",    83: "Bi",    84: "Po",
    85: "At",    86: "Rn",    87: "Fr",    88: "Ra",    89: "Ac",    90: "Th",
    91: "Pa",    92: "U",    93: "Np",    94: "Pu",    95: "Am",    96: "Cm",
    97: "Bk",    98: "Cf",    99: "Es",    100: "Fm",    101: "Md",    102: "No",
    103: "Lr",    104: "Rf",    105: "Db",    106: "Sg",    107: "Bh",    108: "Hs",
    109: "Mt"
}

Element_lmax = {
    "H" : 1, "C": 2, "N": 2, "O" : 2, "F" : 2, "Na": 2, 
    "Mg": 2, "P": 3, "S": 3, "Cl": 3, "Zn": 3,
    "Cu": 3, "Fe": 3, "Ni": 3, "Co": 3
}



#load the input geometry
symbols = []
coords = []
derivative = 0
charge = 0
spin = 1
with open(input_file, 'r') as f:
    line = f.readline()
    arr = line.split()
    natoms = int(arr[0])
    derivative = int(arr[1])
    charge = int(arr[2])
    spin = int(arr[3])

    for i in range(natoms):
        line = f.readline()
        arr = line.split()
        atm = int(arr[0])
        symbol = Element_symbol[atm]
        x = float(arr[1])*0.529177
        y = float(arr[2])*0.529177
        z = float(arr[3])*0.529177
        symbols.append(symbol)
        coords.append([x, y, z])
elements = list(Counter(symbols).keys())
elem_pair = list(itertools.product(elements, elements))


if (derivative == 2):
    print ('The interface is not ready yet for hessian external')
    sys.exit(1)

#check if SK files exist
for elem1, elem2 in elem_pair:
    filename = '{}/{}-{}.{}'.format(SKPATH, elem1, elem2, SKFILE_SUFFIX)
    if (not os.path.exists(filename)):
        print ('SKFILE: {} is not found'.format(filename))
        sys.exit(1)


dirpath = tempfile.mkdtemp()

try:
    with open('{}/dftb.inp'.format(dirpath), 'w') as dftbinp:
        hasForce = False
        hasHessian = False

        #process the template input
        with open(template_input, 'r') as header:
            for line in header:
                if (len(line.strip())==0):
                    break
                if ('FORCE=TRUE' in line):
                    hasForce = True
                    if (derivative == 0):
                        continue
                if ('MD' in line):
                    continue
                if ('OPT' in line):
                    continue
                print(line, end='', file=dftbinp)
            if (derivative > 0 and hasForce == False):
                print('FORCE=TRUE', file=dftbinp)
        print('',file=dftbinp)

        print('Title',file=dftbinp)
        print('',file=dftbinp)

        #process the skfiles
        print(len(elements), file=dftbinp)
        for elem1 in elements:
            print ('{} {}'.format(elem1, Element_lmax[elem1]), file=dftbinp)
            filenames = []
            for elem2 in elements:
                filenames.append('{}/{}-{}.skf'.format(SKPATH, elem1, elem2))
            print(" ".join(filenames), file=dftbinp)
        print('',file=dftbinp)

        #process the geometry
        print ('{} {} {}'.format(natoms, charge, spin), file=dftbinp)
        for i in range(len(symbols)):
            print ('{:2s} {:20.12f}{:20.12f}{:20.12f}'.format(
                symbols[i], *coords[i]
            ), file=dftbinp)
        print('',file=dftbinp)

    with open('{}/dftb.inp'.format(dirpath), 'r') as dftbinp:
        print(''.join(dftbinp.readlines()))

    print ('Running DC-DFTB-K program in {}'.format(dirpath))
    res = subprocess.run([DCDFTBK_PATH], cwd=dirpath)


    #print and parse
    grads = []
    energy = 0.0
    g09out = open(output_file, 'w')
    with open('{}/dftb.out'.format(dirpath), 'r') as dftbout:
        for line in dftbout:
            match = re.match("    Final .* Energy =", line)
            if (match is not None):
                arr = line.split()
                energy = float(arr[4])
                
            if ('Atom           F(x)                F(y)                F(z)' in line):
                grads = []
                next(dftbout)
                for i in range(natoms):
                   line = next(dftbout)
                   arr = line.split()
                   grads.append([-float(arr[1]),-float(arr[2]),-float(arr[3])])
                   
        print('{:20.12f}{:20.12f}{:20.12f}{:20.12f}'.format(
                    energy, 0.0, 0.0, 0.0
                ), file=g09out)
        for line in grads:
            print('{:20.12f}{:20.12f}{:20.12f}'.format(
                    *line 
                ), file=g09out)       
#        print('{:20.12f}{:20.12f}{:20.12f}{:20.12f}'.format(
#                    energy, 0.0, 0.0, 0.0
#                ))
#        for line in grads:
#            print('{:20.12f}{:20.12f}{:20.12f}'.format(
#                    *line 
#                ))       
    g09out.close()         

finally:
    pass
    #shutil.rmtree(dirpath)

