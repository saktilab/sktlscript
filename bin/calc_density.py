#!/usr/bin/env python3

import pymatgen as mg
import sys
from scipy.constants import N_A

mol = mg.Molecule.from_file(sys.argv[1])
mass = sum([x.specie.atomic_mass for x in mol.sites])
mass = mass/N_A

a1 = 0
a2 = 0
a3 = 0 
with open(sys.argv[2], 'r') as f:
	for line in f:
		if "CRYST1" in line:
			arr = line.split()
			a1 = float(arr[1])
			a2 = float(arr[2])
			a3 = float(arr[3])		
V = a1*a2*a3*10**-24 # Volume dalam cm^3
dens = mass/V # Massa jenis dalam g/cm^3

print("Density = {}".format(dens))
