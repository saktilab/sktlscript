#!/usr/bin/env python

import pymatgen as pg
from typing import Union

def loadgen(filename) :#-> Union[pg.Structure,pg.Molecule]:
    elements = []
    coords = []
    with open(filename, 'r') as f:
        line = next(f)
        arr = line.split()
        nat = int(arr[0])
        geomtype = arr[1]
        line = next(f)
        symbols = line.split()
        for _ in range(nat):
            line = next(f)
            arr = line.split()
            elements.append(symbols[int(arr[1])-1])
            coords.append(list(map(float, arr[2:5])))
        if geomtype.upper() == 'C':
            return pg.Molecule(elements, coords=coords)
        else:
            if geomtype.upper() == 'S' or geomtype.upper() == 'F':
                lattice_vectors = []
                line = next(f)
                for _ in range(3):
                    line = next(f)
                    vec = list(map(float, line.split()[:3]))
                    lattice_vectors.append(vec)
                if (geomtype.upper() == 'S'):
                    return pg.Structure(lattice=pg.Lattice(lattice_vectors), species=elements, coords=coords, coords_are_cartesian=True)
                else:
                    return pg.Structure(lattice=pg.Lattice(lattice_vectors), species=elements, coords=coords, coords_are_cartesian=False)
        