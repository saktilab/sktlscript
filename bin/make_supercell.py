#!/usr/bin/env python3

import sys
import argparse
import pymatgen as mg
import random
import numpy as np
import math

"""
This script is for making supercell using DC-DFTB-MD input file

The format of the geometry is the one being used in the DC-DFTB-MD code,
which has the format as follows::

    nat charge spin
    Symbol CoordX CoordY CoordZ
    TV Lat1X Lat1Y Lat1Z
    TV Lat2X Lat2Y Lat2Z
    TV Lat3X Lat3Y Lat3Z
    empty line

Usage::

    make_supercell.py -v 2x2x3 in.dcdftb > out.dcdftb # make the supercell 2x2x3 
    
"""

parser = argparse.ArgumentParser(description='Supercell maker')
parser.add_argument('-v', '--vector', default='1x1x1', type=str, help='supercell', required=True)
parser.add_argument('-c', '--cell', action='store_true')
parser.add_argument('geometry', metavar='geometry', type=str)

opt = parser.parse_args(sys.argv[1:])

factors = list(map(int, opt.vector.split('x')))
if len(factors)!= 3:
    print ('cell definition: AxBxC')
    sys.exit(1)



nat = 0
charge = 0
spin = 1
coords = []
symbols = []
lattice = None

dcdftbmd_input = False
keyword_section = []
title_section = []
parameter_section = []

with open(opt.geometry, 'r') as fin:
    line = next(fin)

    if ('=' in line): #full dcdftbmd input
        dcdftbmd_input = True
        keyword_section.append(line)
        line = next(fin)
        while line.strip():
            keyword_section.append(line)
            line = next(fin)
        line = next(fin)
        while line.strip():
            title_section.append(line)
            line = next(fin)
        line = next(fin)
        while line.strip():
            parameter_section.append(line)
            line = next(fin)

        line = next(fin)
    arr = line.split()
    nat, charge, spin = map(int, arr[0:3])
    for line in fin:
        if line.strip().startswith('TV'):
            if lattice is None:
                lattice = []
            lattice.append( list(map(float, line.split()[1:4])) )
        else:
            arr = line.split()
            if len(arr) < 4:
                continue
            symbols.append(arr[0])
            coords.append(list(map(float, arr[1:4])))
    if len(coords) != nat:
        raise RuntimeError('Number of atom in the geomtry file is not consistent.')
    if (lattice is None or len(lattice) != 3):
        raise RuntimeError('Only 3D PBC (3 TVs) is currently supported.')

mole = mg.Structure(lattice=mg.Lattice(lattice), species=symbols, coords=coords, coords_are_cartesian=True, to_unit_cell=False)
if opt.cell:
        lens, angles = mole.lattice.lengths_and_angles
        new_lens = (lens[0]*factors[0], lens[1]*factors[1], lens[2]*factors[2])
        new_lattice = mg.Lattice.from_lengths_and_angles(new_lens, angles)
        
        
        new_species = []
        new_coords = []
        for x in range(factors[0]):
            for y in range(factors[1]):
                for z in range(factors[2]):
                    for site in mole.sites:
                        new_coord = [site.coords[0]+x*new_lens[0], site.coords[1]+y*new_lens[1], site.coords[2]+z*new_lens[2]]
                        new_species.append(site.specie)
                        new_coords.append(new_coord)
                        
        mole = mg.Structure(new_lattice, new_species, new_coords, coords_are_cartesian=True, to_unit_cell=True)

else:
    mole.make_supercell(factors)

if (dcdftbmd_input):
    print(''.join(keyword_section))
    print(''.join(title_section))
    print(''.join(parameter_section))

print('{} {} {}'.format(len(mole.sites), charge, spin))
for site in mole.sites:
    print('{:<4s} {:13.6f} {:13.6f} {:13.6f}'.format(str(site.specie), *site.coords))
for i in range(3):
    print('{:<4s} {:13.6f} {:13.6f} {:13.6f}'.format('TV', *mole.lattice.matrix[i]))
print("")
