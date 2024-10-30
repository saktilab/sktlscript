#!/usr/bin/env python3

import sys

with open (sys.argv[1], 'r') as f:
    num_atoms = int(next(f))
    charge = []
    atom_serial = []
    next(f)
    atom_count = 0
    for line in f:
        if atom_count < num_atoms:
            arr = line.split()
            atom_serial.append(arr[0])
            charge.append(float(arr[1]))
            atom_count +=1
        
print(atom_serial)


    

    