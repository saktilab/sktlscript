#!/usr/bin/env python3

###Script to modify or adding adsorbate to the surface. 
import pymatgen as mg
import random
import pymatgen.core.surface as sf
import argparse 
import sys

#Remove annoying warnings!
import warnings

parser = argparse.ArgumentParser(description='Build the slab from the bulk structure')
warnings.filterwarnings("ignore")
parser.add_argument('-i', '--input', type=str, required=True, help='Input primary structure of adsorbent (surface in POSCAR format)')
parser.add_argument('-a', '--adsorbate', type=str, required=True, help='Structure of adsorbate (in xyz format)')
parser.add_argument('-d', '--distance', type=float, default=1.0, help='Distance between adsorbent(surface) and adsorbate atoms (in Angstrom)')
# parser.add_argument('-n', '--addnum', type=int, default=1.0, help='Number of adsorbate')
parser.add_argument('-l', '--label', type=str, required=True, help='Label of active sites')
parser.add_argument('-c', '--coverage', type=float,  help='Coverage in the number of atom/Angstroms square')



opt = parser.parse_args(sys.argv[1:])


def add_molecule (struct, mole_ads, ads_index):
    dist = opt.distance
    for ind in ads_index:
        newmole = mole_ads.copy()
        coords = struct.sites[ind].coords
        newmole.translate_sites(vector=[coords[0],coords[1],coords[2]+dist])
        for atom in newmole.sites:
            struct.append(atom.specie, atom.coords, coords_are_cartesian=True)
    
# struct = mg.Structure.from_file(opt.input)
struct = sf.Structure.from_file(opt.input)
mole_ads = mg.Molecule.from_file(opt.adsorbate)
f = open(opt.label, 'r')
for label in f:
    sites_label = [int(x) for x in label.split()]
addnum = int(round(opt.coverage*struct.lattice.a*struct.lattice.b))
print('Number of adsorbates to be added:{}'.format(addnum))
ads = random.sample(sites_label,addnum)
add_molecule(struct,mole_ads,ads)
elements = struct.types_of_specie
element = [str(x) for x in elements]
sorted_struct = struct.get_sorted_structure(key=lambda x: element.index(str(x.specie)))
sorted_struct.to(fmt='poscar', filename='adsorbed_{}.vasp'.format(opt.coverage))