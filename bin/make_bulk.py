#!/usr/bin/env python3
from ase.build import bulk
import sys
import argparse
from ase.io import write,read 
import pymatgen.io.vasp
import numpy as np


def round_up_to_odd(f):
    return int(np.ceil(f)//2*2+1)

parser = argparse.ArgumentParser(description='Build simple bulk structure based on ASE database')
parser.add_argument('-e', '--element', type=str, help='Metal element for creating a bulk structure')
parser.add_argument('-t', '--type', type=str, help='Type of unit cell, e.g. fcc, bcc, sc, hcp, diamond, zincblende, rocksalt, cesiumchloride, fluorite, or wurtzite')
parser.add_argument('-a', '--latt_a', type=float, help='Lattice constant for a-axis')
parser.add_argument('-c', '--latt_c', type=float, help='Lattice constant for c-axis')
parser.add_argument('-covera', '--covera', type=float, help='ratio c/a for hcp unit-cell')
parser.add_argument('-ortho', '--ortho', type=bool, help='Set to true if you want an orthorhombic cell')
parser.add_argument('-cubic', '--cubic', type=bool, help='Set to true if you want a cubic cell')

opt = parser.parse_args(sys.argv[1:])

cell = bulk(opt.element, opt.type, a=opt.latt_a, covera=opt.covera, orthorhombic=opt.ortho, cubic=opt.cubic)
latt_vec = cell.cell

vec_a = [float(i) for i in latt_vec[0]]
vec_b = [float(i) for i in latt_vec[1]]
vec_c = [float(i) for i in latt_vec[2]]    


write('in.gen', cell, format='gen')
write('POSCAR', cell)


# kptsa = round_up_to_odd(50.0/vec_a[0])
# kptsb = round_up_to_odd(50.0/vec_b[1])
# kptsc = round_up_to_odd(50.0/vec_c[2])
# kpoints_obj = pymatgen.io.vasp.Kpoints.gamma_automatic(kpts=(kptsa, kptsb, kptsc))
# kpoints_obj.write_file('KPOINTS')

