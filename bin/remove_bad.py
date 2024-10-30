#!/usr/bin/env python3
import sys
import pymatgen as mg

struct = mg.Structure.from_file(sys.argv[1],format='vasp')
