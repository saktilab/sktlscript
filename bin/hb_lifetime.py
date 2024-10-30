#!/usr/bin/env python3

import MDAnalysis as md 
from MDAnalysis.analysis.waterdynamics import HydrogenBondLifetimes
import pymatgen as mg
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis.analysis.hbonds.hbond_analysis as hb

parser = argparse.ArgumentParser(description='Hydrogen life time calculator')
parser.add_argument('-i', '--input', type=str, help='Input trajectory')
parser.add_argument('-s1', '--selection1',type=str, help='Selection of hydrogen bond acceptor')
parser.add_argument('-s2', '--selection2',type=str, help='Selection of hydrogen bond donor')
parser.add_argument('-l', '--latt', type=str, help='File containing lattice vector file, here you can use dftb.inp of the DC-DFTB-K input')
parser.add_argument('-t0', '--t0',type=int, help='Frame origin')
parser.add_argument('-tf', '--tf', type=int, help='Final frame of origin')
parser.add_argument('-dtmax', '--dtmax', type=int, help='Number of windows created for correlation.')



opt = parser.parse_args(sys.argv[1:])

lattice = []

with open(opt.latt, 'r') as f:
    for line in f:
        if 'TV' in line:
            lattice.append(list(map(float, line.split()[1:4])))

box = mg.Lattice(lattice)
(lx, ly, lz), (a, b, c) = box.lengths_and_angles

u = md.Universe(opt.input, format='xyz')
u.dimensions=np.array([lx, ly, lz, a, b, c], dtype=np.float32)

selection1 = opt.selection1
selection2 = opt.selection2


t0 = opt.t0
tf = opt.tf
dtmax = opt.dtmax

analysis = HydrogenBondLifetimes(u, selection1, selection2, t0, tf, dtmax)

analysis.run()

time = 0
for HBLc, HBLi in analysis.timeseries:
    print(" {time} {HBLc} {time} {HBLi} ".format(time=time, HBLc=HBLc, HBLi=HBLi))
    time += 1

plt.figure(1, figsize=(18,6))

plt.subplot(121)
plt.xlabel('time')
plt.ylabel('HBLc')
plt.title('HBL continuos')
plt.plot(range(0,time), [column[0] for column in analysis.timeseries])

plt.subplot(122)
plt.xlabel('time')
plt.ylabel('HBLi')
plt.title('HBL intermitent')
plt.plot(range(0,time), [column[1] for column in analysis.timeseries])

plt.show()
