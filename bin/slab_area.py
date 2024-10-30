#!/usr/bin/env python3
import pymatgen as mg
import pymatgen.core.surface as sf 
import pymatgen.symmetry.analyzer as sa 
import numpy as np
import argparse 
import sys
import re

#Remove annoying warnings!
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='Build the slab from the bulk structure')
parser.add_argument('POSCAR', metavar='POSCAR', type=str)
parser.add_argument('-hkl', '--miller', type=str, required=True, help='Miller index hkl')
parser.add_argument('-t','--thick', type=float, default=1, help='Thickness of the slab')
parser.add_argument('-v','--vacuum', type=float, default=3, help='Thickness of the vacuum')
parser.add_argument('-c','--center', type=bool, default=True, help='Centering the slab')
parser.add_argument('-u','--unit', type=bool, default=True, help='Define slab in the unit of plane(s)')
parser.add_argument('-p','--primitive', type=bool, default=False, help='Convert to primitive cell')
parser.add_argument('-r','--reorient', type=bool, default=True, help='Reorient the slab orientation')
parser.add_argument('-a', '--super_a', type=int, default=1, help='length of supercell in a axis')
parser.add_argument('-b', '--super_b', type=int, default=1, help='length of supercell in b axis')
parser.add_argument('-n', '--nsearch', type=int, default=None, help='Number of maximum normal search')
###Additional option not included in pymatgen:
parser.add_argument('-s', '--sorting', type=bool, default=False, help='Sorted structure')
parser.add_argument('-sym', '--sga', type=bool, default=False, help='Convert the structure to the standard defined in doi:10.1016/j.commatsci.2010.05.010')
parser.add_argument('-e', '--elements', type=str, required=True, help='Element(s) in the structure, please enter without space, e.g. AlO means that the structure contains Al and O atoms')

opt = parser.parse_args(sys.argv[1:])



bulk = mg.Structure.from_file(opt.POSCAR)

mill = [int(i) for i in str(opt.miller)]


if (opt.sga == True):
    sga = sa.SpacegroupAnalyzer(bulk)
    conv_bulk = sga.get_conventional_standard_structure()
    slab = sf.SlabGenerator(conv_bulk, [mill[0],mill[1],mill[2]], opt.thick, opt.vacuum, center_slab=opt.center, in_unit_planes=opt.unit, primitive=opt.primitive, reorient_lattice=opt.reorient, max_normal_search=opt.nsearch)
else:
    slab = sf.SlabGenerator(bulk, [mill[0],mill[1],mill[2]], opt.thick, opt.vacuum, center_slab=opt.center, in_unit_planes=opt.unit, primitive=opt.primitive, reorient_lattice=opt.reorient, max_normal_search=opt.nsearch)




bulk = mg.Structure.from_file(opt.POSCAR)

mill = [int(i) for i in str(opt.miller)]


if (opt.sga == True):
    sga = sa.SpacegroupAnalyzer(bulk)
    conv_bulk = sga.get_conventional_standard_structure()
    slab = sf.SlabGenerator(conv_bulk, [mill[0],mill[1],mill[2]], opt.thick, opt.vacuum, center_slab=opt.center, in_unit_planes=opt.unit, primitive=opt.primitive, reorient_lattice=opt.reorient)
else:
    slab = sf.SlabGenerator(bulk, [mill[0],mill[1],mill[2]], opt.thick, opt.vacuum, center_slab=opt.center, in_unit_planes=opt.unit, primitive=opt.primitive, reorient_lattice=opt.reorient)

slab_coord = slab.get_slab()
slab_coord.make_supercell([[opt.super_a,0,0],[0,opt.super_b,0],[0,0,1]])

print('# Surface area(in meter square):',slab_coord.surface_area*1e-20)
print('# Surface area(in Angstrom square):',slab_coord.surface_area)
