#!/usr/bin/env python3

import os
import numpy as np 
import sys
from shutil import move
import argparse

parser = argparse.ArgumentParser(description='Prepare data set for various c/a ratio of hcp structures')
parser.add_argument('-e', '--element', type=str, help='Metal element for creating a bulk structure')
parser.add_argument('-t', '--type', type=str, help='Type of unit cell, e.g. fcc, bcc, sc, hcp, diamond, zincblende, rocksalt, cesiumchloride, fluorite, or wurtzite')
parser.add_argument('-a', '--latt_a', type=float, help='Lattice constant for a-axis')
parser.add_argument('-start', '--start', type=float, help='Starting point of c/a ratios')
parser.add_argument('-end', '--end', type=float, help='End point of c/a ratios')
parser.add_argument('-step', '--step', type=float, help='Step of c/a ratios')

opt = parser.parse_args(sys.argv[1:])


coveras = np.arange(float(opt.start), float(opt.end), float(opt.step))
coveras = [float(i) for i in coveras]

for covera in coveras:
    covera = str(round(covera,3))
    if not os.path.exists(covera):
        os.makedirs(covera)
    os.chdir(covera)
    element = str(opt.element)
    unitcell = str(opt.type)
    lattice = str(opt.latt_a)
    cmd = '~/bin/make_bulk.py -e ' + element + ' -t ' + unitcell + ' -a ' + lattice +  ' -ortho True ' + ' -covera ' + covera
    os.system(cmd)
    os.chdir("..")
        
    
