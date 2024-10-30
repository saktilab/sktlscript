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
bonds_bohr = []
for ifile in input_files:
    with open(ifile, 'r') as f:
        struct = pymatgen.Molecule.from_str(f.read(), fmt='xyz')
        for ind, site1 in enumerate(struct.sites):
            if (str(site1.specie) != elem1):
                continue
            for ind2, site2 in enumerate(struct.sites):
                # if (ind <= ind2):
                if (str(site2.specie) != elem2):
                    continue
                r = site1.distance(site2)
                if (r > cutoff):
                    continue
                bonds.append(r)

print('Size:    ', len(bonds))
print('Min:     {:20.12f} {:20.12f}'.format( min(bonds), min(bonds)/0.529177))
print('Max:     {:20.12f} {:20.12f}'.format( max(bonds), max(bonds)/0.529177))
print('Mean:    {:20.12f} {:20.12f}'.format( statistics.mean(bonds), statistics.mean(bonds)/0.529177))
print('Median:  {:20.12f} {:20.12f}'.format( statistics.median(bonds), statistics.median(bonds)/0.529177))
print('Std.Dev: {:20.12f} {:20.12f}'.format( statistics.stdev(bonds), statistics.stdev(bonds)/0.529177))

bonds_bohr = np.array(bonds)/0.529177
plot = plt.hist(bonds_bohr, bins=100)
plt.xlabel("Diatomic Distance [a.u.]")
plt.ylabel("Counts")
plt.title('{}-{}'.format(elem1, elem2))
# plt.show()
plt.savefig('{}-{}.eps'.format(elem1, elem2), format="eps")
