#!/usr/bin/env python3

import sys
import numpy as np
from scipy.constants import k

cv = []
height = []

with open(sys.argv[1], 'r') as f:
    for line in f:
        if ("Coordinate" in line):
            arr = line.split()
            cv.append(float(arr[2]))
        if ("Gaussian height" in line):
            arr = line.split()
            # height.append(float(arr[3])*1059.70/298*float(sys.argv[2])) # in kT unit         
            height.append(float(arr[3])) # in hartree unit
            
cv = np.array(cv, dtype='float')
height = np.array(height, dtype='float')
colvar_inp = np.column_stack((cv,height))

np.savetxt("COLVAR", colvar_inp)