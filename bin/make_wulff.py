#!/usr/bin/env python3
from ase.cluster import wulff_construction
from ase.io import read, write
import sys
import argparse

import warnings

parser = argparse.ArgumentParser(description='Build cluster based on Wulff theorem\n not that at present this only work for a cubic lattice')
warnings.filterwarnings("ignore")
parser.add_argument('-s1', '--s1', type=int, default=100, help='Miller index for surface 1')
parser.add_argument('-s2', '--s2', type=int, default=110, help='Miller index for surface 2')
parser.add_argument('-s3', '--s3', type=str, default=111, help='Miller index for surface 3')
parser.add_argument('-lc', '--lc', type=float, help='The lattice constant.')
parser.add_argument('-natom','--natom', type=int, default=100, help='Approximate number of atoms')
parser.add_argument('-e','--element', type=str, help='Element')
parser.add_argument('-u','--unitcell', type=str, help='Type of unit cell, e.g. fcc, bcc, or sc')
parser.add_argument('-r','--rounding', type=str, default="closest", help='Specifies what should be done if no Wulff construction corresponds to exactly the requested number of atoms. Should be a string, either “above”, “below” or “closest” (the default), meaning that the nearest cluster above or below - or the closest one - is created instead')
parser.add_argument('-fmt','--format', type=str, default='xyz', help='format of the output file')
parser.add_argument('-es1','--es1', type=float,  help='Surface energy for surface 1')
parser.add_argument('-es2','--es2', type=float,  help='Surface energy for surface 2')
parser.add_argument('-es3','--es3', type=float,  help='Surface energy for surface 3')
parser.add_argument('-i','--input', type=str,  help='Input file contains miller index and surface energy in two columns')


opt = parser.parse_args(sys.argv[1:])

# mill1 = [int(i) for i in str(opt.s1)]
# mill2 = [int(i) for i in str(opt.s2)]
# mill3 = [int(i) for i in str(opt.s3)]

f = open(opt.input, 'r')
esurf = []
miller = []
for line in f:
    miller.append(line[0])
    esurf.append(line[1])

# surfaces = [(mill1[0],mill1[1],mill1[2]),(mill2[0],mill2[1],mill2[2]), (mill3[0],mill3[1],mill3[2])]
# esurf = [opt.es1, opt.es2, opt.es3]
#Initialization of lattice constant
lc = opt.lc
size = opt.natom
print(miller)
print(esurf)
atoms = wulff_construction(opt.element, miller, esurf, size, opt.unitcell, opt.rounding, latticeconstant=lc)

if (opt.format == 'xyz'):
    write('wulff.xyz', atoms)
else:
    write('wulff.vasp', atoms)
