#!/usr/bin/env python3


import pymatgen as mg
import sys


struct = mg.Structure.from_file(sys.argv[1])
print (struct.lattice.volume/(0.529177**3.0)/float(len(struct.sites)))

