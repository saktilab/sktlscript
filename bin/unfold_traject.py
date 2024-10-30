#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import pymatgen 
from pymatgen.util import coord

parser = argparse.ArgumentParser(description='Compute MSD')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')
parser.add_argument('-s', '--start', type=int, default=0, help='StartTime')
parser.add_argument('-g', '--groups', type=str, required=True, help='Molecules, use name1*nat1*nmole1:name2*nat2*nmole2')
opt = parser.parse_args(sys.argv[1:])

index_groups = []

entries = opt.groups.split(':')
for entry in entries:
    arr = entry.split('*')
    if (len(arr) != 3):
        print('Format error for groups')
        parser.print_help()
        sys.exit(1)
    index_groups.append( (arr[0], int(arr[1]), int(arr[2])) )
# print(index_groups)


traj_file = opt.traject
# print('Starting step: ', opt.start)

symbols = []
symbol_list = []
nat_list = []
lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = pymatgen.Lattice(lattice)

#read the first geometry from the trajectory file
with open(traj_file, 'r') as f:
    nat = int(f.readline())
    info = f.readline()
    arr = info.split()
    t1 = float(arr[3])
    for i in range(0, nat):
        line = f.readline()
        arr = line.split()
        symbols.append(arr[0])
        if (len(symbol_list) == 0):
            symbol_list.append(arr[0])
            nat_list.append(1)
        else:
            if (arr[0] == symbol_list[-1]):
                nat_list[-1] += 1
            else:
                symbol_list.append(arr[0])
                nat_list.append(1)
    nat = int(f.readline())
    info = f.readline()
    arr = info.split()
    t2 = float(arr[3])
    dt = t2-t1

groups_nat = sum([k[1]*k[2] for k in index_groups])
if (groups_nat != nat):
    raise ValueError('Wrong groups spec. The # of atoms are not consistent')

# print('Groups:')
# for k, nat, nmole in index_groups:
#     print ('    Name: {}, # of moles: {}, # of atoms per mole: {}'.format(k, nat, nmole))
# print("")
    


def load_xyz(lattice, f):
    nat = int(next(f))
    title = next(f)
    species = []
    coords = []
    for i in range(nat):
        line = next(f)
        arr = line.split()
        species.append(arr[0])
        coords.append(list(map(float, arr[1:4])))
    mole = pymatgen.Structure(lattice=lattice, species=species, coords=coords, coords_are_cartesian=True)
    return title, mole

last_coords = []
with open(traj_file, 'r') as f:
    step = 0
    while(True):
        try:
            step += 1
            
            title, structure = load_xyz(lattice, f)
            
            print(len(structure.sites))
            print(title.rstrip())

            atom_index = 0
            group_index = 0
            cur_coords = []
            for group in index_groups:
                # print('Compute COM for group {}'.format(group))
                for mole_index in range(group[2]):
                    # print('Compute COM for molecule {}'.format(mole_index))
                    local_indexs = [atom_index + x for x in range(group[1])]
                    # print(local_indexs)
                    
                    atom1 = structure[local_indexs[0]]

                    if (step == 1):   
                        cur_coords.append(atom1.coords)
                    
                        for at2 in local_indexs[1:]:
                            atom2 = structure[at2]
                            dis_, image_ = atom1.distance_and_image(atom2)
                            cur_coords.append(structure.lattice.get_cartesian_coords(atom2.frac_coords+image_))                       
                        
                    else:
                        for at in local_indexs:
                            atom = structure[at]
                            vec = pymatgen.util.coord.pbc_shortest_vectors(structure.lattice, 
                                structure.lattice.get_fractional_coords(last_coords[at]),
                                atom.frac_coords
                                )
                            cur_coords.append(last_coords[at]+vec[0][0]) 
                        
                    atom_index = local_indexs[-1]+1
            
            last_coords = cur_coords

            for ind, item in enumerate(cur_coords):
                print('{:<5s} {:20.12f} {:20.12f} {:20.12f}'.format(str(structure.sites[ind].specie), *item))
                
            
            


        except StopIteration as e:
            break

    



