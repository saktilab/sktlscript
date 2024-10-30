#!/usr/bin/env python

import sys
import pymatgen as mg

logfile = sys.argv[1]

nat = 0
step = 1

TV = []
symbols = []
coords = []


with open(logfile, 'r') as f:
    for line in f:
        if 'Total number of atoms' in line:
            arr = line.split()
            nat = int(arr[5])
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
        
struct = mg.Structure(mg.Lattice(TV), species=symbols, coords=coords, coords_are_cartesian=True)
print(struct.to(fmt='poscar'))