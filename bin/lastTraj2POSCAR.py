#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Grep the last geometry from XYZ trjatory \ncombining with the lattice')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')

opt = parser.parse_args(sys.argv[1:])

traj_file = opt.traject
lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(line.split()[1:4])

if (len(lattice)!=3):
    raise ValueError('Number of lattice vectors is not 3.')

lines = []
with open(traj_file, 'r') as f:
    for line in f:
        lines = []
        nat = int(line)
        lines.append(line)
        title = next(f)
        lines.append(title)
        for i in range(nat):
            line = next(f)
            lines.append(line)

li = iter(lines)
nat = int(next(li))
title=next(li)
symbols = []
symbol_list = []
nat_list = []
coords = []
for i in range(nat):
    line = next(li)
    arr = line.split()
    symbols.append(arr[0])
    coords.append(arr[1:4])
    if (len(symbol_list) == 0):
        symbol_list.append(arr[0])
        nat_list.append(1)
    else:
        if (symbol_list[-1] == arr[0]):
            nat_list[-1] += 1
        else:
            symbol_list.append(arr[0])
            nat_list.append(1)

print(title.strip())
print(1.0)
for line in lattice:
    print(' '.join(line))
print(' '.join(symbol_list))
print(' '.join(map(str, nat_list)))
print('Cart')
for coord, sym in zip(coords, symbols):
    print('{:20.12f} {:20.12f} {:20.12f}   {:3s}'.format(*list(map(float,coord)), sym))


