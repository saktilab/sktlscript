#!/usr/bin/env python3
import pymatgen as mg
import pymatgen.core.surface as sf 
import pymatgen.analysis as pa
import pymatgen.symmetry.analyzer as sa 
import numpy as np
import argparse 
import sys
import re

#Remove annoying warnings!
import warnings
warnings.filterwarnings("ignore")


structure = mg.Structure.from_file(sys.argv[1])
area = structure.lattice.a * structure.lattice.b

print('# Surface area(in meter square):',area*1e-20)
print('# Surface area(in Angstrom square):',area)
