#!/usr/bin/env python3
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.core.bonds import get_bond_length, obtain_all_bond_lengths
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Compute Coordination Number of an Element in a Cluster System')
parser.add_argument('-i', '--input', type=str, required=True, help='Trajectory input file')
parser.add_argument('-s', '--site', type=int, required=True, help='Index of element')
parser.add_argument('-dt', '--time_step', type=float, default=1.0, help='Time step of MD simulation')
parser.add_argument('-p', '--print_freq', type=int, default=10, help='Printing frequency in MD simulation')
opt = parser.parse_args(sys.argv[1:])

print('#Time [ps]', '#Coordination Number')
with open(opt.input, 'r') as f:
    step = 0
    for line in f:
        # coords = []
        x = []
        y = []
        z = []
        nat = int(line)
        info = next(f)
        arr = info.split()
        for i in range(0, nat):
            line = next(f)
            arr = line.split()
            # coords.append(np.array(list(map(float, arr[1:4])),dtype=np.double))
            x.append(np.array(arr[1], dtype=np.double))
            y.append(np.array(arr[2], dtype=np.double))
            z.append(np.array(arr[3], dtype=np.double))
        active_sites = [int(opt.site)]
        CN = 0
        bulk_dist = 2.703
        step+=1
        for i in active_sites:
            CN = 0
            for j in range(0, nat):
                if (i == j):
                    continue
                dist = ((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)**0.5
                if (dist <= 3.3):
                    CN = CN + (2*(1-(dist/bulk_dist)**6)/(1-(dist/bulk_dist)**12))
            print(step*opt.print_freq*opt.time_step/1000, CN)
        