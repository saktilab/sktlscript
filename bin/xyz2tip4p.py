#!/usr/bin/env python3

import pymatgen as mg
import sys
import numpy as np

inputfile = sys.argv[1]

mole = mg.Molecule.from_file(inputfile)

nat = len(mole.sites)
n_waters = int(nat/3)

final_specie = []
final_coords = []

#assuming the first line is O and the following 2 lines are H
for ind in range(n_waters):
    i = ind*3
    mid_of_h = (mole.sites[i+1].coords + mole.sites[i+2].coords)/2.0
    vec = (mole.sites[i].coords - mid_of_h)
    vec = (vec / np.linalg.norm(vec))*0.15
    final_coord = mole.sites[i].coords+vec
    final_specie.append('OH2') #OH2
    final_specie.append('H1') #H1
    final_specie.append('H2') #2
    final_specie.append('OM') #OM
    final_coords.append(mole.sites[i].coords)
    final_coords.append(mole.sites[i+1].coords)
    final_coords.append(mole.sites[i+2].coords)
    final_coords.append(final_coord)

#xyz output
#print(len(final_specie))
#print(' ')
#for i in range(len(final_coords)):
#    print('{} {:20.12f} {:20.12f} {:20.12f}'.format(final_specie[i], *final_coords[i]))

### pdb output
for i in range(n_waters):
    seq_id = i+1
    for j in range(4):
        atom_id = i*4+j+1

        print ('ATOM  {:5d} {:>4s} TIP4Q{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      W1'
        .format( atom_id,
            final_specie[atom_id-1],
            seq_id,
            *final_coords[atom_id-1],
            1.0, 0.0)
        )
