#!/usr/bin/env python3
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.core.bonds import get_bond_length, obtain_all_bond_lengths
import sys

with open(sys.argv[1], 'r') as f:
    F = f.readlines()
    Fn = int(F[0])
    Fm = Molecule.from_file(sys.argv[1])
    active_site = [int(sys.argv[2])] 
    S = 0
    i = 0
    CN = 0
    bulk_dist = 2.703
    # bulk_dist = 4
    for na in active_site:
        CN = 0
        for nb in range(0, Fn):
            if (na == nb):
                continue
            Fd = Molecule.get_distance(Fm, na, nb)
            if (Fd <= 3.3):
                CN = CN + (2*(1-(Fd/bulk_dist)**6)/(1-(Fd/bulk_dist)**12))
        print(CN)