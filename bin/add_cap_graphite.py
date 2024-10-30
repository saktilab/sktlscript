#!/usr/bin/env python3
import sys
import numpy as np

symbol = []
x = []
y = []
z = []
with open (sys.argv[1], "r") as f:
    Natom = next(f)
    next(f)
    for line in f:
        arr = line.split()
        symbol.append(arr[0])
        x.append(float(arr[1]))
        y.append(float(arr[2]))
        z.append(float(arr[3]))
        
coord = np.array(list(zip(x,y,z)))
coord = np.sort(coord, axis=0)

ala = np.nonzero(coord == -0.6957500)
print(ala)