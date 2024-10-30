#!/usr/bin/env python3
import pymatgen as mg
import pymatgen.analysis.structure_analyzer as sa
import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np


struct1 = mg.Structure.from_file(sys.argv[1])
struct2 = mg.Structure.from_file(sys.argv[2])

compare = sa.RelaxationAnalyzer(struct1,struct2)

distance = compare.get_percentage_bond_dist_changes() 

atom1 = []
atom2 = []
dist_percentage = []

for key, values in distance.items():
    atom1.append(key)
    for key, subvalues in values.items():
        atom2.append(key)
        dist_percentage.append(subvalues)

MAE = np.mean(np.abs(dist_percentage))
# MSE = np.mean(dist_percentage)
print("Bond distance difference [%]: {}".format(MAE*100))

    
    
