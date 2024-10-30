#!/usr/bin/env python3

import sys
import argparse
from pymatgen.core import Molecule, Structure
import random
import numpy as np
import math

"""
This script is for adding/removing hydrogen(s) to/from oxygens to make
hydroniums/hydroxides.

The format of the geometry is the one being used in the DC-DFTB-MD code,
which has the format as follows::

    nat charge spin
    Symbol CoordX CoordY CoordZ
    TV Lat1X Lat1Y Lat1Z
    TV Lat2X Lat2Y Lat2Z
    TV Lat3X Lat3Y Lat3Z
    empty line

The lattice vectors (TV) are optional if no PBC applies.

Usage::

    add_hydrogen.py -n 1 in.dcdftb > out.dcdftb # add one hydrogen to the system
    add_hydrogen.py -n -2 in.dcdftb > out.dcdftb # remove 2 hydrogens from the system

Note:
Only tested on PBC systems.

Powered by Chien-Pin Chou

version 20190313: fix a bug while providing only the geometry part of the input
version 20190312-2: initial


All Copyright Reserved

"""

def eulerAnglesToRotationMatrix(theta) :

    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])



    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])

    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])


    R = np.dot(R_z, np.dot( R_y, R_x ))

    return R



def rotation_matrix(angle, u):
    sina = math.sin(angle)
    cosa = math.cos(angle)

    R = np.zeros((3,3))
    R[0][0] = cosa + u[0]*u[0]*(1.0-cosa)
    R[0][1] = u[0]*u[1]*(1.0-cosa)-u[2]*sina
    R[0][2] = u[0]*u[2]*(1.0-cosa)+u[1]*sina

    R[1][0] = u[0]*u[1]*(1.0-cosa)+u[2]*sina
    R[1][1] = cosa + u[1]*u[1]*(1.0-cosa)
    R[1][2] = u[2]*u[1]*(1.0-cosa)-u[0]*sina

    R[2][0] = u[2]*u[0]*(1.0-cosa)-u[1]*sina
    R[2][1] = u[1]*u[2]*(1.0-cosa)+u[0]*sina
    R[2][2] = cosa + u[2]*u[2]*(1.0-cosa)

    return R




parser = argparse.ArgumentParser(description='Hydrogen Adder')
parser.add_argument('-n', '--numHydrogen', default=1, type=int, help='Number of hydrogen(s) to be added/removed')
parser.add_argument('-d', '--debug', action='store_true', help='Write the result to debug.vasp/debug.xyz')
parser.add_argument('geometry', metavar='geometry', type=str)

opt = parser.parse_args(sys.argv[1:])

if (opt.numHydrogen == 0):
    print('Nothing to do.')
    sys.exit(0)

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
    if (lattice is not None and len(lattice) != 3):
        raise RuntimeError('Only 3D PBC (3 TVs) is currently supported.')

if lattice is None:
    mole = Molecule(species=symbols, coords=coords)
else:
    mole = Structure(lattice=mg.Lattice(lattice), species=symbols, coords=coords, coords_are_cartesian=True)

neutral_water_indexes = []
hydroxide_indexes = []
hydronium_indexes = []
all_oxygen_indexes = []

for ind, site in enumerate(mole.sites):
    if str(site.specie) != 'O':
        continue
    neighbors = mole.get_neighbors(site, 2.0)
    neighbors_h = [x for x in neighbors if (str(x[0].specie) == 'Al' or str(x[0].specie) == 'Si') and x[1] < 1.8]
    if len(neighbors_h) == 4:
        neutral_water_indexes.append(ind)
    else:
        if len(neighbors_h) == 1 and len(neighbors) == 1: #hydroxide
            hydroxide_indexes.append(ind)
        elif len(neighbors_h) == 3 and len(neighbors) == 3:
            hydronium_indexes.append(ind)

random.shuffle(neutral_water_indexes)
random.shuffle(hydroxide_indexes)
random.shuffle(hydronium_indexes)

if opt.numHydrogen > 0 :
    all_oxygen_indexes = hydroxide_indexes + neutral_water_indexes
else:
    all_oxygen_indexes = hydronium_indexes + neutral_water_indexes


if (abs(opt.numHydrogen) > len(all_oxygen_indexes)):
    raise ValueError('Number of hydrogens to add/remove is more than the number of oxygens')

oxygen_to_alter = all_oxygen_indexes[:abs(opt.numHydrogen)]


final_symbols = []
final_coords = []

if opt.numHydrogen > 0:
    #add hydrogen
    processed_hydrogens = 0
    new_h_coords = []
    for ind in oxygen_to_alter:

            neighbors = [x for x in mole.get_neighbors(mole.sites[ind], 1.2, include_index=True) if str(x[0].specie) == 'H' ]
            if len(neighbors) == 2:
                # very simple algorithm that produces a "flat" H3O+
                h1_coord = neighbors[0][0].coords
                h2_coord = neighbors[1][0].coords
                o_coord = mole.sites[ind].coords
                v1 = h1_coord-o_coord
                v2 = h2_coord-o_coord
                new_h_coord = o_coord-0.8*(v1+v2)
                new_h_coords.append(new_h_coord.tolist())
                processed_hydrogens += 1
            elif len(neighbors) == 1:
                 # very simple algorithm that produces a "flat" H3O+
                h1_coord = neighbors[0][0].coords
                o_coord = mole.sites[ind].coords
                v1 = h1_coord-o_coord

                theta = [random.random()*math.pi for i in range(3)]
                rot_mat = eulerAnglesToRotationMatrix(theta)
                v2 = np.dot(rot_mat, v1)
                v_n = np.cross(v1, v2)
                v_n = v_n/np.linalg.norm(v_n)

                rot_mat = rotation_matrix(72/180.0*math.pi, v_n)

                new_vec = np.dot(rot_mat, v1)
                new_h_coord = o_coord-new_vec
                new_h_coords.append(new_h_coord.tolist())
                processed_hydrogens += 1
            else:
                continue


    if (processed_hydrogens != len(oxygen_to_alter)):
        raise RuntimeError('Some of the selected oxygen does not belong to neutral water.')

    final_symbols = symbols + ['H']*processed_hydrogens
    final_coords = coords + new_h_coords

else:
    #remove hydrogen
    processed_hydrogens = 0
    removed_h_indexes = []
    for ind in oxygen_to_alter:
            neighbors = [x for x in mole.get_neighbors(mole.sites[ind], 1.2, include_index=True) if str(x[0].specie) == 'H' ]

            if len(neighbors) == 1:
                continue

            remove_index = random.choice([x[2] for x in neighbors])
            removed_h_indexes.append(remove_index)
            processed_hydrogens += 1


    if (processed_hydrogens != len(oxygen_to_alter)):
        raise RuntimeError('Some of the selected oxygen does not belong to neutral water.')

    final_symbols = symbols
    final_coords = coords
    for ind in sorted(removed_h_indexes, reverse=True):
        del final_symbols[ind]
        del final_coords[ind]

if (opt.numHydrogen < 0):
    processed_hydrogens = -processed_hydrogens

if (opt.debug):
    if isinstance(mole, mg.core.Structure):
        new_mole = mg.Structure(lattice=mole.lattice, species=final_symbols, coords=final_coords, coords_are_cartesian=True)
        new_mole.sort()
        new_mole.to(filename='debug.vasp', fmt='poscar')
    else:
        new_mole = mg.Molecule(species=final_symbols, coords=final_coords)
        new_mole.to(filename='debug.xyz', fmt='xyz')


if (dcdftbmd_input):
    print(''.join(keyword_section))
    print(''.join(title_section))
    print(''.join(parameter_section))

print('{} {} {}'.format(len(final_coords), charge+processed_hydrogens, spin))
for sym, xyz in zip(final_symbols, final_coords):
    print('{:<4s} {:13.6f} {:13.6f} {:13.6f}'.format(sym, *xyz))
if lattice is not None:
    for i in range(3):
        print('{:<4s} {:13.6f} {:13.6f} {:13.6f}'.format('TV', *lattice[i]))
print("")
