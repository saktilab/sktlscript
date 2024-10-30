#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import numpy.linalg as npa
import pymatgen as mg
import copy
import pymatgen.util.coord as cu

parser = argparse.ArgumentParser(description='Compute Grotthus shuttling')
parser.add_argument('input', metavar='input', type=str)
parser.add_argument('-l', '--latt', type=str, required=True)
parser.add_argument('-s', '--init_step', type=int, default=-1)
opt = parser.parse_args(sys.argv[1:])
fout = open("gro.dat", "w")
lattice = []
with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)

def get_vector(box, fcoord1, fcoord2):
    dv = cu.pbc_shortest_vectors(box, fcoord2, fcoord1)[0][0]
    return dv



input_file = opt.input


num_pros = 0
max_charge_atoms = []

# previous_coords = None

# print ('#{:>9s} {:>8s} {:>8s} {:>10s}'.format('Time(ps)', '# Grotthus', 'Index', 'Distance'))
print ('#{:>9s} {:>8s} {:>8s}'.format('Time(ps)', '# Grotthus', 'Index'),file=fout)
with open(input_file, 'r') as f:

    interrupted = False

    previous_previous_index = 0
    previous_index = 0
    previous_backward_index = 0
    step = opt.init_step -1
    for line in f:
        step += 1
        if (line.startswith('#')):
            continue
        arr = line.split()
        # step = int(arr[0])
        time = float(arr[0])
        index = int(arr[2])
        # coords = list(map(float, arr[4:]))
        
        dis = 0.0
        if (index == -1):
            interrupted = True
            # print ('{:10.4f} {:>8d} {:>8d} {:10.4f}'.format(time, num_pros, index, dis))
            print ('{:10.4f} {:>8d} {:>8d}'.format(time, num_pros, index),file=fout)
            continue
        
        if (step == 0):
            previous_previous_index = index
            previous_index = index
            previous_backward_index = index
            # previous_coords = coords.copy()
        elif (step == 1):
            if (index != previous_index):
                num_pros += 1
                previous_index = index
                # dv = get_vector(box, previous_coords, coords)            
                # dis = npa.norm(dv)
                # previous_coords = coords.copy()
        else:
            if (index == previous_index):
                pass
            else:
                # dv = get_vector(box, previous_coords, coords)    
                # dis = npa.norm(dv)
                # print ('dis', index, previous_index, dis)
                # previous_coords = coords.copy()
                if (index != previous_index and index == previous_previous_index):
                    if (previous_index == previous_backward_index):
                        num_pros += 1
                    else:
                        num_pros += -1
                        previous_backward_index = index
                else:
                    num_pros += 1
                previous_previous_index = previous_index
                previous_index = index

            # print ('{:10.4f} {:>8d} {:>8d} '.format(time, num_pros, max_charge_atoms[step]))

        print ('{:10.4f} {:>8d} {:>8d}'.format(time, num_pros, index),file=fout)

