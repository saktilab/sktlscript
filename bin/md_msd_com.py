#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import pymatgen as mg
import pymatgen.util.coord as cu

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
print('Starting step: ', opt.start)

symbols = []
symbol_list = []
nat_list = []
lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)

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

print('Groups:')
for k, nat, nmole in index_groups:
    print ('    Name: {}, # of atoms per mole: {}, # of moles: {}'.format(k, nat, nmole))
print("")
    

first = True
first_coords = []


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
    mole = mg.Structure(lattice=lattice, species=species, coords=coords, coords_are_cartesian=True)
    return mole


fout = open('msd.out', 'w')
format_str = '{:>15s}(A^2)'*len(index_groups)

print ('#Time(ps) ', format_str.format(*[k for k, nat, nmole in index_groups]), file=fout)

first_coms = None
last_coms = None
image_lasts = None

times = []
msds = []

with open(traj_file, 'r') as f:
    step = 0
    while(True):
        try:
            step += 1
            
            structure = load_xyz(lattice, f)
            atom_index = 0
            group_index = 0
            com_step = []
            for group in index_groups:
                # print('Compute COM for group {}'.format(group))
                for mole_index in range(group[2]):
                    # print('Compute COM for molecule {}'.format(mole_index))
                    local_indexs = [atom_index + x for x in range(group[1])]
                    # print(local_indexs)
                    if (len(local_indexs) == 1):
                        com_step.append(structure.sites[local_indexs[0]].frac_coords)
                    else:
                        atom1 = structure[local_indexs[0]]
                        com_local = atom1.specie.data['Atomic mass']*atom1.frac_coords
                        total_mass = atom1.specie.data['Atomic mass']
                        for at2 in local_indexs[1:]:
                            atom2 = structure[at2]
                            dis_, image_ = atom1.distance_and_image(atom2)
                            com_local += atom2.specie.data['Atomic mass']*(atom2.frac_coords+image_)
                            total_mass += atom2.specie.data['Atomic mass']
                        com_local /= total_mass
                        com_step.append(com_local)
                    atom_index = local_indexs[-1]+1
            
            if (step == 1):
                first_coms = com_step.copy()
                last_coms = com_step.copy()
                image_lasts = [(0,0,0)]*len(last_coms)
            else:
                msd = []
                image_cur_list = []
                
                for first_com, last_com, cur_com, image_last in zip(first_coms, last_coms, com_step, image_lasts):
                    image_current = structure.lattice.get_distance_and_image(last_com, cur_com)[1]
                    image_final = image_last + image_current
                    final_dis = structure.lattice.get_distance_and_image(first_com, cur_com, jimage=image_final)[0]
                    msd.append(final_dis*final_dis)
                    image_cur_list.append(image_final)
                
                image_lasts = image_cur_list
                last_coms = com_step.copy()

                msd_group = []
                for i in range(len(index_groups)):
                    items = msd[0:index_groups[i][2]]
                    # print(len(items))
                    msd = msd[index_groups[i][2]:]
                    msd_group.append(sum(items)/len(items))
                format_str = '{:<10.5f} ' + '{:20.12f}'*len(index_groups)
                print (format_str.format(step*dt/1000.0, *msd_group), file=fout)
                
            
            print('Step {} done'.format(step))


        except StopIteration as e:
            break

    



