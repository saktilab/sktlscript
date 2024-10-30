#!/usr/bin/env python3
import pymatgen as mg
import pymatgen.core.surface as sf 
import pymatgen.symmetry.analyzer as sa 
import numpy as np
import argparse 
import sys
import os
import pymatgen.io.vasp as pv


#Remove annoying warnings!
import warnings

parser = argparse.ArgumentParser(description='###Build the slab from the bulk structure###')

warnings.filterwarnings("ignore")
parser.add_argument('POSCAR', metavar='POSCAR', type=str)
parser.add_argument('-hkl', '--miller', type=str, required=True, help='Miller index hkl')
parser.add_argument('-t','--thick', type=float, default=1, help='Thickness of the slab')
parser.add_argument('-v','--vacuum', type=float, default=3, help='Thickness of the vacuum')
parser.add_argument('-c','--center', type=str, default=True, help='Centering the slab')
parser.add_argument('-u','--unit', type=bool, default=True, help='Define slab in the unit of plane(s)')
parser.add_argument('-p','--primitive', type=str, default=True, help='Convert to primitive cell')
parser.add_argument('-r','--reorient', type=str, default=True, help='Reorient the slab orientation')
parser.add_argument('-a', '--super_a', type=int, default=1, help='length of supercell in a axis')
parser.add_argument('-b', '--super_b', type=int, default=1, help='length of supercell in b axis')
parser.add_argument('-n', '--nsearch', type=int, default=None, help='Number of maximum normal search')
###Additional option not included in pymatgen:
parser.add_argument('-s', '--sorting', type=bool, default=True, help='Sorted structure')
parser.add_argument('-sym', '--sga', type=bool, default=False, help='Convert the structure to the standard defined in doi:10.1016/j.commatsci.2010.05.010')
parser.add_argument('-reduce', '--reduce', type=bool, default=False, help='Whether or not the slabs will be orthogonalized')
parser.add_argument('-shift', '--shift', type=float, default=0, help='Shift the c-direction of the slab')

opt = parser.parse_args(sys.argv[1:])

def round_up_to_odd(f):
    return int(np.ceil(f) // 2 * 2 + 1)



bulk = mg.core.Structure.from_file(opt.POSCAR)

mill = [int(i) for i in str(opt.miller)]


if (opt.sga == True):
    sga = sa.SpacegroupAnalyzer(bulk)
    conv_bulk = sga.get_conventional_standard_structure()
    slab = sf.SlabGenerator(conv_bulk, [mill[0],mill[1],mill[2]], opt.thick, opt.vacuum, center_slab=opt.center, in_unit_planes=opt.unit, primitive=opt.primitive, reorient_lattice=opt.reorient, max_normal_search=opt.nsearch, lll_reduce=opt.reduce)
else:
    slab = sf.SlabGenerator(bulk, [mill[0],mill[1],mill[2]], opt.thick, opt.vacuum, center_slab=opt.center, in_unit_planes=opt.unit, primitive=opt.primitive, reorient_lattice=opt.reorient, max_normal_search=opt.nsearch, lll_reduce=opt.reduce)

slab_coord = slab.get_slab(shift=opt.shift)
slab_coord.make_supercell([[opt.super_a,0,0],[0,opt.super_b,0],[0,0,1]])
f = open("area.dat", "w")
elements = slab_coord.types_of_specie
element = [str(x) for x in elements]

f.write("Surface area in angstrom: {}".format(slab_coord.surface_area))
if (opt.sorting == True): 
    slabsort=sf.Slab.get_sorted_structure(slab_coord,key=lambda x: element.index(str(x.specie)))
    slabsort.to(filename='slab.{}.vasp'.format(int(opt.thick)),fmt='poscar')
else:
    slab_coord.to(filename='slab.{}.vasp'.format(int(opt.thick)),fmt='poscar')

kptsa = round_up_to_odd(50.0/slab_coord.lattice.a)
kptsb = round_up_to_odd(50.0/slab_coord.lattice.b)

kpoints = pv.Kpoints.gamma_automatic(kpts=(kptsa,kptsb,1))
kpoints.write_file("KPOINTS")
