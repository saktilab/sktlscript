#!/usr/bin/env python3

import sys
import pymatgen as mg

with open(sys.argv[1], 'r') as f:
    input_str = mg.Structure.from_str(f.read(), fmt='poscar')

    print('{} 0 1'.format(len(input_str.sites)))
    for site in input_str.sites:
        print('{:<3s} {:20.12f} {:20.12f} {:20.12f}'.format(str(site.specie), *site.coords))
    for i in range(3):
        print('TV  {:20.12f} {:20.12f} {:20.12f}'.format(*input_str.lattice.matrix[i]))
    print() 
