#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Compute Proton transfer')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')
parser.add_argument('-s', '--steps', type=int, nargs='+', help='time steps')
parser.add_argument('-if', '--index-file', type=str, help='Index file to load and analyzed')
opt = parser.parse_args(sys.argv[1:])

traj_file = opt.traject

index2print = []
with open(opt.index_file, 'r') as f:
    for line in f:
        line_list = list(map(int, line.split()))
        index2print += line_list

lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

step = 0
ind = 0
with open(traj_file, 'r') as f:
    for line in f:
        coords = []
        nat = int(line)
        info = next(f)
    
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            if (i in index2print):
                coords.append( np.array(list(map(float, arr[1:4])),dtype=np.double))
            
        if (step in opt.steps):
            with open('{}.xyz'.format(ind+1), 'w') as fout:
                print(len(coords), file=fout)
                print(lattice[0][0], file=fout)
                print(lattice[1][1], file=fout)
                print(lattice[2][2], file=fout)
                for line in coords:
                    print('{} {} {}'.format(*line), file=fout)
            ind += 1
        step += 1
        
        