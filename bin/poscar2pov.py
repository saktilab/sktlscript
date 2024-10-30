#!/usr/bin/env python3
from ase.io.vasp import read_vasp
from ase.io import read, write

import sys
import argparse
parser = argparse.ArgumentParser(description='Convert xyz,POSCAR, or gen files into povray input.')
parser.add_argument('-i', '--input', type=str, help='input file geometry')
parser.add_argument('-r', '--rotate', type=str, help='rotation term to rotate the resulting figure. Example of rotation format: 1x 1y 1z.')

opt = parser.parse_args(sys.argv[1:])

atoms = read_vasp(opt.input)

write('structure.pov', atoms, rotation=opt.rotate)    
