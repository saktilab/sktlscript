#!/usr/bin/env python3
from ase.io import read, write
import os
import sys
import argparse
parser = argparse.ArgumentParser(description='Convert xyz files into povray input.')
parser.add_argument('-i', '--input', type=str, help='input file geometry')
parser.add_argument('-r', '--rotate', type=str, help='rotation term to rotate the resulting figure. Example of rotation format: 1x 1y 1z.')

opt = parser.parse_args(sys.argv[1:])
os.system("echo {} > geom.smi".format(opt.input))
os.system("obabel geom.smi -O geom.xyz --gen3D")

atoms = read("geom.xyz")

write('structure.pov', atoms, rotation=opt.rotate)    
