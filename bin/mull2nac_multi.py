#!/usr/bin/env python3

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='DC-DFTB-K mulliken to NAC')
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('mulliken', metavar='mulliken', type=str)
opt = parser.parse_args(sys.argv[1:])

#number of shells
n_orbs = {
    'Pt': 3,
    'Br': 3,
    'O': 2,
    'C': 2,
    'N': 2,
    'H': 1,
    'S': 3,
    'Li': 2,
    'F': 2,
    'P': 3
}

val = {
    'Li': 1,
    'H': 1,
    'C': 4,
    'O': 6,
    'N': 5,
    'S': 6,
    'F': 7,
    'P': 5
}
traj_file = opt.traject
mull_file = opt.mulliken


symbols = []
symbol_list = []
nat_list = []


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


#parse the mulliken charge file and compute the net atomic charge(NAC)
step = 1
#avg_charge = []
fout = open("nac.dat", "w") 
with open(mull_file, 'r') as f:
    for line in f:
        header = line.split()
        # total_orbs = header[0]
        nat = header[0] 
        stepline = next(f)
        print(nat, file=fout)
        print(stepline, end='', file=fout)
        arr = stepline.split()
        step += 1

        for ind, elem in enumerate(symbols):
            norb = n_orbs[elem]
            nac = 0.0
            for j in range(0, norb):
                line = next(f)
                arr = line.split()                
                nac += float(arr[1])
            nac = val[elem] - nac  
            print('{:6d} {:20.12f}'.format(ind+1, nac), file=fout)
fout.close()
