#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import Rbf
import matplotlib

x = []
y = []
z = []
with open(sys.argv[1], "r") as f:
    for line in f:
        x.append(float(line))
with open(sys.argv[2], "r") as f:
    for line in f:
        y.append(float(line))
with open(sys.argv[3], "r") as f:
    for line in f:
        z.append(float(line))

x = np.array(x, dtype=np.double)
y = np.array(y, dtype=np.double)
z = np.array(z, dtype=np.double)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi,yi)

rbf = Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

# plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto', cmap='winter_r')
fig = plt.figure()
contour = plt.contour(xi, yi, zi, 20, colors='k')
plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
contour_filled = plt.contourf(xi, yi, zi, 20)
plt.colorbar(contour_filled)
plt.xlabel("N-O distance [Angstroms]",fontsize=20)
plt.ylabel("Rh-N-O angle [Deg.]",fontsize=20)
# plt.scatter(x, y, c=z)
plt.tick_params(direction='in')
# plt.colorbar()
plt.show()

