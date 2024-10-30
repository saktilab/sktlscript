#!/usr/bin/env python3
# Script for calculating binding energy in kcal/mol unit
import sys

complex = str(sys.argv[1])
protein = str(sys.argv[2])
ligand = str(sys.argv[3])


with open(complex, 'r') as f:
    for line in f:
        if "Final" in line:
            arr = line.split()
            E_complex = float(arr[4])

with open(protein, 'r') as f:
    for line in f:
        if "Final" in line:
            arr = line.split()
            E_protein = float(arr[4])

with open(ligand, 'r') as f:
    for line in f:
        if "Final" in line:
            arr = line.split()
            E_ligand = float(arr[4])

E_bind = (E_complex - E_ligand - E_protein)*627.509

print("Binding energy = {:.2f} kcal/mol".format(E_bind))