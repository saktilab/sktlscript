#!/usr/bin/env python3


import pymatgen as mg
import sys
import xml.etree.ElementTree as ET
import os


path = sys.argv[1]

struct = mg.Structure.from_file(os.path.join(path, 'CONTCAR'))
volume = struct.lattice.volume/(0.529177**3.0)/float(len(struct.sites))

doc = ET.ElementTree()
doc.parse(os.path.join(path,'vasprun.xml'))
root = doc.getroot()
energy = float(root.find('calculation/energy/i').text)/float(len(struct.sites))/27.2112

print ('{:20.12f} {:20.12f}'.format(volume, energy))
