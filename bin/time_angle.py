#!/usr/bin/env python3
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.core.bonds import get_bond_length, obtain_all_bond_lengths
import sys
import numpy as np

with open(sys.argv[1], 'r') as f:
    step = 0
    for line in f:
        coords = []
        # x = []
        # y = []
        # z = []
        nat = int(line)
        info = next(f)
        arr = info.split()
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            coords.append(np.array(list(map(float, arr[1:4])),dtype=np.double))
            # x.append(np.array(arr[1], dtype=np.double))
            # y.append(np.array(arr[2], dtype=np.double))
            # z.append(np.array(arr[3], dtype=np.double))
        a = int(sys.argv[2])
        b = int(sys.argv[3])
        c = int(sys.argv[4])
        ba = coords[a]-coords[b]
        bc = coords[c]-coords[b]

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba)*np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        print (np.degrees(angle))
        