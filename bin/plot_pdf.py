#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

x = []
y = []

with open(sys.argv[1], "r") as f:
    for line in f:
        arr = line.split()
        x.append(float(arr[1]))
with open(sys.argv[2], "r") as f:
    for line in f:
        arr = line.split()
        y.append(float(arr[1]))

x = np.array(x)
y = np.array(y)

p = x[None,:]*y[:,None]

plt.figure(dpi=300)

plt.imshow(p,cmap='gnuplot2', origin='lower',extent=[2,3,0,1],interpolation='None')
plt.colorbar()
plt.show()
