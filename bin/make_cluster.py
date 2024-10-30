#!/usr/bin/env python3
import ase
# from ase.cluster.icosahedron import Icosahedron
import ase.cluster as clus
from ase.io import read, write
import sys
import argparse 

import warnings
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(description='Build cluster structures based on element information')
parser.add_argument('-e', '--element', type=str, help='Element' )
parser.add_argument('-t', '--type', type=str, help='Type of cluster, e.g. octahedron, decahedron, and  icosahedron' )
parser.add_argument('-s', '--shell', type=int, help='Number of shell of the cluster' )
parser.add_argument('-lc', '--lc', type=float, help='Lattice constant of the bulk structure' )
opt = parser.parse_args(sys.argv[1:])

if (opt.type == "icosahedron"):
    struct = clus.Icosahedron(opt.element,opt.shell, opt.lc)
    struct.write("{}_{}_{}.xyz".format(opt.element,opt.type,opt.shell))
if (opt.type == "cubic"):
    struct = clus.cubic(opt.element,opt.shell, opt.lc)
    struct.write("{}_{}_{}.xyz".format(opt.element,opt.type,opt.shell))
if (opt.type == "octahedron"):
    struct = clus.Octahedron(opt.element,opt.shell, opt.lc)
    struct.write("{}_{}_{}.xyz".format(opt.element,opt.type,opt.shell))
if (opt.type == "decahedron"):
    struct = clus.Decahedron(opt.element,opt.shell, opt.lc)
    struct.write("{}_{}_{}.xyz".format(opt.element,opt.type,opt.shell))
if (opt.type == "tetrahedron"):
    struct = clus.Tetrahedron(opt.element,opt.shell, opt.lc)
    struct.write("{}_{}_{}.xyz".format(opt.element,opt.type,opt.shell))

# elif (opt.type == "decahedron"):
#     struct = clus.Decahedron(opt.element,opt.shell, opt.lc)
#     struct.write("{}_{}_{}.vasp".format(opt.element,opt.type,opt.shell))
# write(struct, 'Pt147.vasp', format='vasp')

