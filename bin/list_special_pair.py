#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pymatgen as mg

parser = argparse.ArgumentParser(description='List the distance')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-l', '--lattice', metavar='lattice', type=str)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--index', nargs='+', type=int, help='Index to be analyzed')
group.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')

opt = parser.parse_args(sys.argv[1:])

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

traj_file = opt.traject
latt_file = opt.lattice

lattice_vec = []
with open(latt_file, 'r') as f:
    for line in f:
        if ('TV' in line):
            lattice_vec.append( list(map(float,line.split()[1:4])) )

Latt = mg.Lattice(lattice_vec)

list_pairs = []


if (opt.index):
    if (len(opt.index) % 2 != 0):
        raise ValueError('The length of index should be even.')
    list_pairs = list(pairwise(opt.index))
elif (opt.index_file):
    with open(opt.index_file, 'r') as f:
        for line in f:
            line_list = list(map(int, line.split()[0:2]))
            list_pairs.append(line_list)

pbc = False
if (len(lattice_vec) == 3):
    pbc = True


symbols = []

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
    nat = int(f.readline())
    info = f.readline()
    arr = info.split()
    t2 = float(arr[3])
    dt = (t2-t1)/1000.0

final_pairs = []

print ()
print ('Pairs to be printed:')
for li in list_pairs:
    
    if (li[0] < len(symbols) and li[1] < len(symbols)):
        print ('{:5d} {:5d}'.format(*li))
        final_pairs.append([li[0], li[1]])
print ()


fout = open('list_spec_pair.dat', 'w')

step = 0
# format_str = ' {:>15s}'*len(final_pairs)
# headers = ['{}({})-{}({})'.format(symbols[i], i, symbols[j], j) for i, j in final_pairs ]
# print(('{:>10s} ' + format_str).format('#Time(ps)', *headers), file=fout)
with open(traj_file, 'r') as f:
    for line in f:
        time = dt*step
        # print("Time: {} ps".format(time), file=fout)
        step+=1
        nat = int(line)
        next(f)
        coords = []
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            coords.append(list(map(float, list(arr[1:4]))))
        if (pbc):
            my_str = mg.Structure(Latt, species=symbols, coords=coords, coords_are_cartesian=True)
        else:
            my_str = mg.Molecule(species=symbols, coords=coords)
        distances = []
        for ind_1, ind_2 in final_pairs:
            distance = my_str.get_distance(ind_1, ind_2)
            if (distance < 1.2):
                print("{} {} {}".format(time, ind_2, distance), file=fout)
        # format_str = ' {:15.4f}'*len(distances)
        
        # print(('{:>10.4f} ' + format_str + ' {:15.4f}').format(time, *distances, sum(distances)/len(distances)),file=fout)
        # print(('{:>10.4f} ' + format_str).format(time, *distances),file=fout)    
    