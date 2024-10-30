#!/usr/bin/env python3

import sys
import pymatgen.core.structure as mg

struct = mg.Structure.from_file(sys.argv[1])

vac = [int(x)-1 for x in sys.argv[2:]]
struct.remove_sites(vac)

struct.to(fmt='poscar', filename="vacant.vasp")