#!/usr/bin/env python3
import numpy as np
import sys
from matplotlib import pyplot as plt

x = []
y = []
with open(sys.argv[1], "r") as f:
    for line in f:
        arr = line.split()
        x.append(float(arr[0]))
        y.append(float(arr[1]))
for x1, y1 in zip(x,y):
    plt.plot(x1,y1,'ro')

z = np.polyfit(x, y, 4)
f = np.poly1d(z)

for x1 in np.linspace(0, 0.6, 100):
    plt.plot(x1, f(x1), 'b+')
plt.axis([0,1,0,3])
plt.show()