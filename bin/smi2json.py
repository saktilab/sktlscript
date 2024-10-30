#!/usr/bin/env python3
import pybel
import json

molecule = pybel.readstring("smi", "O=C1C2=C(N=CN2C)N(C(=O)N1C)C")
molecule.make3D()

with open("caffeine.json", "w") as out_file:
    out_file.write(molecule_to_json(molecule))
