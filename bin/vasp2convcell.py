#!/usr/bin/env python3

import sys
import pymatgen as mg

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer 

filename = sys.argv[1] 

with open(filename, 'r') as f:    
    struct = mg.Structure.from_str(f.read(), fmt='poscar')
    sa = SpacegroupAnalyzer(struct)
    conv_cell = sa.get_conventional_standard_structure()
    print(conv_cell.to(fmt='poscar'))
