#!/usr/bin/env python

import sys
import pymatgen as mg

logfile = sys.argv[1]

readTV = False
if (len(sys.argv)==3):
    readTV = True



nat = 0
step = 1

TV = []
symbols = []
coords = []
spin = 1
charge = 0

if (readTV):
    with open(sys.argv[2], 'r') as f:
        for line in f:
            if 'TV' in line:
                TV.append(list(map(float, line.split()[1:4])))


with open(logfile, 'r') as f:
    for line in f:
        if 'Total number of atoms' in line:
            arr = line.split()
            nat = int(arr[5])
        if ('Spin multiplicity') in line:
            arr = line.split()
            spin = int(arr[3])
        if ('Charge of system') in line:
            arr = line.split()
            charge = int(arr[4])
        if 'Molecular coordinate (angs)' in line:
            next(f)
            next(f)
            next(f)
            next(f)
            
            step += 1
            symbols = []
            coords = []
            for i in range(nat):
                line = next(f)
                arr = line.split()
                symbols.append(arr[0])
                coords.append(list(map(float, arr[1:4])))
            

        if 'Lattice vector (angs)' in line:
            next(f)
            next(f)
            next(f)
            next(f)
            
            TV = []
            for i in range(3):
                line = next(f)
                arr = line.split()
                TV.append(list(map(float, arr[1:4])))
        
# struct = mg.Structure(mg.Lattice(TV), species=symbols, coords=coords, coords_are_cartesian=True)
# print(struct.to(fmt='poscar'))

print('{} {} {}'.format(nat, charge, spin))
for symb, coord in zip(symbols, coords):
    print('{:3s} {:20.12f} {:20.12f} {:20.12f}'.format(symb, *coord))

for lat in TV:
    print('{:3s} {:20.12f} {:20.12f} {:20.12f}'.format('TV', *lat))

print()
