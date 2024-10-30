#!/usr/bin/env python3

import sys
import numpy as np
from scipy.constants import k

cv = []
height = []

with open(sys.argv[1], 'r') as f:
    for line in f:
        if ("Gaussian height" in line):
            arr = line.split()
            # height.append(float(arr[3])*1059.70/298*float(sys.argv[2])) # in kT unit         
            height.append(float(arr[3])) # in hartree unit


with open(sys.argv[2], 'r') as f:
    for line in f:
        arr = line.split()
        cv.append(arr[1])



cv = np.array(cv, dtype='float')
height = np.array(height, dtype='float')
# Modifying Gaussian height so it has the same dimension as CV
mod_height = []
for i in range(int(int(sys.argv[3])/int(sys.argv[4]))):
    mod_height.append(0)
for i in height:
    mod_height.append(i)
    for j in range(int(int(sys.argv[3])/int(sys.argv[4])-1)):
        mod_height.append(0)
mod_height = np.array(mod_height, dtype='float')
for i in range(int(int(sys.argv[3])/int(sys.argv[4])-1)):
    mod_height = np.delete(mod_height,-1)

for i in mod_height:
    print(i)
colvar_inp = np.column_stack((cv,mod_height))

np.savetxt("COLVAR", colvar_inp)