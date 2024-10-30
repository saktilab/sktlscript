#!/usr/bin/env python3
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)
import sys

structure = parser.get_structure("pdb", sys.argv[1])

# atom = structure[0]["X"][100]["CA"]
atoms = structure.get_atoms()
for atom in atoms:
    print(atom)
