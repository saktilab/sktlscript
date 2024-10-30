#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

data = []
with open(sys.argv[1], "r") as f:
    for line in f:
        arr = line.split()
        data.append(float(arr[0]))

data = np.array(data)

fig, ax = plt.subplots()
ax.tick_params(axis="y",direction="in")
ax.tick_params(axis="x",direction="in")


n, bins, patches = plt.hist(x=data, bins='auto', color='green', alpha=0.7, rwidth=0.85,edgecolor='black')
# plt.grid( alpha=0.75)
plt.xlabel('Dispersion Contribution [kcal/mol]')
plt.ylabel('# Trajectory')
maxfreq = n.max()
# plt.ylim(top=np.ceil(maxfreq/10)*10 if maxfreq % 10 else maxfreq + 10)

plt.savefig("histo.pdf", format='pdf')