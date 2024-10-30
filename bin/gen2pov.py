#!/usr/bin/env python3
from ase.io import read, write
import pymatgen as mg

#Load the structure
structure = mg.Structure.from_file('POSCAR')
