#!/usr/bin/env python3

from ase.cluster.octahedron import Octahedron
from ase.cluster.icosahedron import Icosahedron
from ase.cluster.decahedron import Decahedron

from pymatgen.symmetry.analyzer import PointGroupAnalyzer

import numpy as np
import pymatgen as mg
from pymatgen.core import Structure, Molecule, Lattice
import os.path
import sys

def make_decahedron(Symbol, param, lc=None):
    p, q, r = param
    atoms = Decahedron(Symbol, p=p, q=q, r=r, latticeconstant=lc)
    return Molecule([Symbol]*len(atoms.positions),atoms.positions)

def make_icosahedron(Symbol, param, lc=None):
    p = param
    atoms = Icosahedron(Symbol, noshells=p, latticeconstant=lc)
    return Molecule([Symbol]*len(atoms.positions),atoms.positions)

def make_octahedron(Symbol, param, lc=None):
    l = param
    atoms = Octahedron(Symbol, length=l, latticeconstant=lc)
    return Molecule([Symbol]*len(atoms.positions),atoms.positions)

def make_trancated_octahedron(Symbol, param, length=None, lc=None):
    c = param
    l = c*3+1
    if (length is not None ): # and length >= c*3+1):
        l = length
    atoms = Octahedron(Symbol, length=l, cutoff=c, latticeconstant=lc)
    return Molecule([Symbol]*len(atoms.positions),atoms.positions)

def make_cuboctahedron(Symbol, param, lc=None):
    c = param
    l = c*2+1
    atoms = Octahedron(Symbol, length=l, cutoff=c, latticeconstant=lc)
    return Molecule([Symbol]*len(atoms.positions),atoms.positions)

def get_box_size(mole, separation=10.0, interval=5.0):
    diameter = max(map(max, mole.distance_matrix))
    cur_size = 0.0
    while(cur_size < (diameter+separation)):
        cur_size += interval
    lattice = Lattice.from_lengths_and_angles((cur_size, cur_size, cur_size), (90,90,90))
    return lattice

def make_pbc_structure(lattice, mole):
    return Structure(lattice=lattice, species=mole.species, coords=mole.cart_coords, coords_are_cartesian=True)


para_decahedron_indexes = [ (x,y,z) for x in range(1,6) for y in range(1,6) for z in range(0,6) ]
para_icosahedron_indexes = [2, 3, 4, 5, 6, 7, 8, 9]
para_cuboctahedron_indexes = [1, 2, 3, 4, 5, 6, 7, 8, 9]
para_t_octahedron_indexes = [1, 2, 3, 4, 5, 6]
para_octahedron_indexes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

para_icosahedron_indexes = list(range(2, 20))
para_cuboctahedron_indexes = list(range(1,20))

all_types = [
    ('icosa', make_icosahedron, para_icosahedron_indexes),
    ('cubo',  make_cuboctahedron, para_cuboctahedron_indexes),
    ('deca',  make_decahedron, para_decahedron_indexes),
    ('tranocta', make_trancated_octahedron, para_t_octahedron_indexes),
    ('octa', make_octahedron, para_octahedron_indexes),            
]

#all_types = [ 
 #   ('cubo',  make_cuboctahedron, para_cuboctahedron_indexes), 
 #   ('icosa', make_icosahedron, para_icosahedron_indexes),
#]


if __name__ == "__main__":
    Symbol = str(sys.argv[1])
    lattice = float(sys.argv[2])
    pwd = os.curdir
    for name, func, indexes in all_types:
        flist = open('{}.list'.format(name), 'w')
        os.makedirs(name, exist_ok=True)
        os.chdir(name)
        for index in indexes:
            mole = func(Symbol, index, lc=lattice)
#            if (len(mole.sites) > 3000):
#                continue
            print('{}{}'.format(Symbol, len(mole.sites)), file=flist)
            print (name, index, len(mole.sites))
            filename = '{}_{}.xyz'.format(Symbol, len(mole.sites))
            mole.to(fmt='xyz', filename=filename)
                # struct = make_pbc_structure(lattice, mole)
                # filename = '{}{}_{:.3f}.vasp'.format(Symbol, len(struct.sites), lat)
                # print(filename)           
                # struct.to(fmt="poscar", filename=filename)
        os.chdir('../')
        flist.close()
    
