#!/usr/bin/env python3

import sys
import pymatgen
import statistics
import numpy as np
import matplotlib.pyplot as plt


elem1 = sys.argv[1]
elem2 = sys.argv[2]
cutoff = float(sys.argv[3])

input_files = sys.argv[4:]

pairs = [(elem1, elem2), (elem2, elem1)]

bonds = []
for ifile in input_files:
    with open(ifile, 'r') as f:
        struct = pymatgen.Structure.from_str(f.read(), fmt='poscar')
        for ind, site in enumerate(struct.sites):
            neighbors = struct.get_neighbors(site, r=cutoff, include_index=True)
            for neighbor, r, n_ind in neighbors:
                if (n_ind < ind):
                    continue
                if (str(site.specie), str(neighbor.specie)) in pairs:
                    bonds.append(r)

print('Size:    ', len(bonds))
print('Min:     {:20.12f} {:20.12f}'.format( min(bonds), min(bonds)/0.529177))
print('Max:     {:20.12f} {:20.12f}'.format( max(bonds), max(bonds)/0.529177))
print('Mean:    {:20.12f} {:20.12f}'.format( statistics.mean(bonds), statistics.mean(bonds)/0.529177))
print('Median:  {:20.12f} {:20.12f}'.format( statistics.median(bonds), statistics.median(bonds)/0.529177))
print('Std.Dev: {:20.12f} {:20.12f}'.format( statistics.stdev(bonds), statistics.stdev(bonds)/0.529177))

bonds_bohr = np.array(bonds)/0.529177
plot = plt.hist(bonds_bohr, bins=100, density=True)
plt.xlabel("Diatomic Distance [a.u.]")
plt.ylabel("Density")
plt.title('{}-{}'.format(elem1, elem2))
# plt.show()
plt.savefig('{}-{}.eps'.format(elem1, elem2), format="eps")
