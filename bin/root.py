#!/usr/bin/env python3
import pymatgen as mg
import numpy as np
import sys

p = [float(x) for x in sys.argv[1:]]

roots=np.roots(p)
print("The roots are:")
index = 0
for root in roots:
	index += 1
	print("root[{}]: {:4.4f}".format(index, root))


