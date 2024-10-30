#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')
parser.add_argument('-s', '--start', type=int, default=0, help='StartTime')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')
opt = parser.parse_args(sys.argv[1:])


traj_file = opt.traject

index2print = []
if (opt.index_file):
    with open(opt.index_file, 'r') as f:
        for line in f:
            line_list = list(map(int, line.split()))
            index2print += line_list

lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

boxlen = []
for i in range(3):
    boxlen.append(lattice[i][i])

step = 0
with open(traj_file, 'r') as f:
    for line in f:
        coords = []
        nat = int(line)
        info = next(f)
        arr = info.split()
    
        for i in range(0, nat):
            line = next(f)
            if step >= opt.start:
                arr = line.split()
                if (i in index2print):
                    coords.append( np.array(list(map(float, arr[1:4])),dtype=np.double))
            
        if step == opt.start:
            print(len(coords))
            for i in range(3):
                print(boxlen[i])
            for coord in coords:
                print('{} {} {}'.format(*coord))
        
            sys.exit(0)
        step += 1
        
        
