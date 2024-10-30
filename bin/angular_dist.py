#!/usr/bin/env python3

import MDAnalysis as md 
from MDAnalysis.analysis.waterdynamics import AngularDistribution as AD
import pymatgen as mg
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Angular distribution of water molecules')
parser.add_argument('-i', '--input', type=str, help='Input trajectory')
parser.add_argument('-s', '--selection',type=str, help='Selection of water molecules in the trajectory')
parser.add_argument('-b', '--bins', type=int, help='Number of bins')
parser.add_argument('-l', '--latt', type=str, help='File containing lattice vector file, here you can use dftb.inp of the DC-DFTB-K input')

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

analysis = AD(u, opt.selection, opt.bins)
analysis.run()



#Plot the results

plt.figure(1, figsize=(18,6))

plt.subplot(131)
plt.xlabel('cos theta')
plt.ylabel('P(cos theta)')
plt.title('PDF cos theta of OH')
plt.plot([float(column.split()[0]) for column in analysis.graph[0][:-1]], [float(column.split()[1]) for column in analysis.graph[0][:-1]])

plt.subplot(132)
plt.xlabel('cos theta')
plt.ylabel('P(cos theta)')
plt.title('PDF cos theta of HH')
plt.plot([float(column.split()[0]) for column in analysis.graph[1][:-1]], [float(column.split()[1]) for column in analysis.graph[1][:-1]])

plt.subplot(133)
plt.xlabel('cos theta')
plt.ylabel('P(cos theta)')
plt.title('PDF cos theta of dipole')
plt.plot([float(column.split()[0]) for column in analysis.graph[2][:-1]], [float(column.split()[1]) for column in analysis.graph[2][:-1]])

plt.show()