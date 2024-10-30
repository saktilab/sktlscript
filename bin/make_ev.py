#!/usr/bin/env python3

import os
import numpy as np 
import sys
from shutil import move
import argparse

parser = argparse.ArgumentParser(description='Prepare data set for E-V curve with the variation of unit cell')
parser.add_argument('-e', '--element', type=str, help='Metal element for creating a bulk structure')
parser.add_argument('-t', '--type', type=str, help='Type of unit cell, e.g. fcc, bcc, sc, hcp, diamond, zincblende, rocksalt, cesiumchloride, fluorite, or wurtzite')
parser.add_argument('-a', '--latt_a', type=float, help='Lattice constant for a-axis')
parser.add_argument('-i', '--incr', type=float, help='Increment of increasing/decreasing of lattice vector (in percent)')
parser.add_argument('-p', '--npoints', type=int, help='Number of points for constructing the set')
# parser.add_argument('-start', '--start', type=float, help='Starting point of c/a ratios')
# parser.add_argument('-end', '--end', type=float, help='End point of c/a ratios')
# parser.add_argument('-step', '--step', type=float, help='Step of c/a ratios')

opt = parser.parse_args(sys.argv[1:])

start = float(opt.latt_a)-float(opt.incr)*0.5*float(opt.npoints)/100
end = float(opt.latt_a) + float(opt.incr)*0.5*float(opt.npoints)/100
step = float(opt.incr)/100
lattices = np.arange(start, end, step)
lattices = [float(i) for i in lattices]

for lattice in lattices:
    lattice = str(round(lattice,3))
    if not os.path.exists(lattice):
        os.makedirs(lattice)
    os.chdir(lattice)
    element = str(opt.element)
    unitcell = str(opt.type)
    cmd = '~/bin/make_bulk.py -e ' + element + ' -t ' + unitcell + ' -a ' + lattice +  ' -cubic True ' 
    os.system(cmd)
    os.chdir("..")

    
