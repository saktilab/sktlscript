#!/usr/bin/env python3
from ase.cluster import wulff_construction
from ase.io import read, write
from pymatgen.core.surface import Structure, Lattice
from pymatgen.analysis.wulff import WulffShape
import sys
import argparse

import warnings

parser = argparse.ArgumentParser(description='Build cluster based on Wulff theorem\n not that at present this only work for a cubic lattice\n please specifically prepare a file, i.e. "surf_ener.dat" that contains dictionary of (h,k,l): surface energy in J/m^2')
warnings.filterwarnings("ignore")


parser.add_argument('-lc', '--lc', type=float, help='The lattice constant.')
parser.add_argument('-natom','--natom', type=int, default=100, help='Approximate number of atoms')
parser.add_argument('-e','--element', type=str, help='Element')
parser.add_argument('-u','--unitcell', type=str, help='Type of unit cell, e.g. fcc, bcc, or sc')
parser.add_argument('-r','--rounding', type=str, default="closest", help='Specifies what should be done if no Wulff construction corresponds to exactly the requested number of atoms. Should be a string, either “above”, “below” or “closest” (the default), meaning that the nearest cluster above or below - or the closest one - is created instead')
parser.add_argument('-fmt','--format', type=str, default='xyz', help='format of the output file')
parser.add_argument('-i','--input', type=str,  help='Two columns table containing surface energy')
parser.add_argument('-show','--show', type=bool,  help='Showing the detail of area fractions')

opt = parser.parse_args(sys.argv[1:])

data = []
with open(opt.input, 'r') as inf:
        data = eval(inf.read())

surfaces = list(data.keys())
esurf = list(data.values())

f = open(opt.input, 'r')



lc = opt.lc
size = opt.natom

atoms = wulff_construction(opt.element, surfaces, esurf, size, opt.unitcell, opt.rounding, latticeconstant=lc)

if (opt.format == 'xyz'):
    write('wulff.xyz', atoms)
elif (opt.format == 'gen'):
    write('in.gen', atoms) 
else:
    write('wulff.vasp', atoms)

write('wulff.pov', atoms, rotation='10z,-80x')

lattice = Lattice.cubic(lc)

wulffshape = WulffShape(lattice, surfaces, esurf)
fraction = wulffshape.area_fraction_dict
print("#Miller index","#Area fraction")
for i in fraction:
    print(i, fraction[i])
print("shape factor: %.3f, anisotropy: \
%.3f, weighted surface energy: %.3f J/m^2" %(wulffshape.shape_factor,                    wulffshape.anisotropy, wulffshape.weighted_surface_energy))

if(opt.show == True):
    ala = wulffshape.get_plot(color_set='PuBu', grid_off=True, axis_off=True, show_area=True, alpha=1, off_color='red', direction=(1,0,0), bar_pos=(0.75,0.15,0.05,0.65), bar_on=False, legend_on=True, aspect_ratio=(8,8))

    ala.savefig("wulff_area.pdf", format='pdf', orientation = 'landscape', papertype=None, transparent=True, bbox_inches='tight')
   