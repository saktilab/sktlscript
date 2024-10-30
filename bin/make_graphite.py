#!/usr/bin/env python3
# from ase.build import graphene_nanoribbon
from ase.io import read, write
import sys
import argparse

import warnings

# parser = argparse.ArgumentParser(description='Build graphite structure')
# warnings.filterwarnings("ignore")
# parser.add_argument("-w", type=float, help="Width of the graphene layer")
# parser.add_argument("-l", type=float, help="Length of the graphene layer")
# parser.add_argument("-type", type=str, default='armchair', help='Type of graphene layer, i.e. armchair or zigzag default: %(default)s')
# parser.add_argument("-vacuum", type=float, default=3.5, help="Vacuum layer thickness in non periodic direction [Angstroms], default: %(default)s")
# parser.add_argument("-r_cc", type=float, default=1.09, help="Length of C-C bond, default:%(default)s")
# parser.add_argument("-r_ch",type=float, default=1.42, help="Length of C-H bond, default=%(default)s")

# opt = parser.parse_args(sys.argv[1:])

from ase.lattice.hexagonal import *
index1=9
index2=9
alat = 2.45
clat = 3.5
gra = Graphene(symbol = 'C',latticeconstant={'a':alat,'c':clat},size=(index1,index2,5))
write("ala.xyz", gra, format='xyz')
