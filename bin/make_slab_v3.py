#!/usr/bin/env python3
from ase.io import read, write
from ase.build import surface
import argparse
import sys


parser = argparse.ArgumentParser(description='Make a slab structure')
parser.add_argument('-hkl', '--hkl', type=str, required=True, help='Indeks Miller')
parser.add_argument('-n', '--layer', type=int, required=True, help='Number of slab layer')
parser.add_argument('-s', '--structure', type=str, required=True, help='Bulk structure in POSCAR format')
parser.add_argument('-v', '--vacuum', type=int, default=20, help='Number of slab layer')
parser.add_argument('-a', '--a', type=int, required=True, help='a axis dimension')
parser.add_argument('-b', '--b', type=int, required=True, help='b axis dimension')
opt = parser.parse_args(sys.argv[1:])

hkl = [int(x) for x in str(opt.hkl)]

bulk = read(opt.structure)

slab = surface(bulk, (hkl[0],hkl[1],hkl[2]), opt.layer, vacuum=opt.vacuum)


write('slab.vasp', slab*(opt.a,opt.b,1))
write('slab.xyz', slab*(opt.a,opt.b,1))


