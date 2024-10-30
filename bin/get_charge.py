#!/usr/bin/env python3

import sys
import numpy as np
import os


charge=0
with open("{}.mol2".format(sys.argv[1]), 'r') as f:
	next(f)
	next(f)
	next(f)
	next(f)
	next(f)
	next(f)
	next(f)
	next(f)
	# head = next(f).split()
	# natom = int(head[0])
	for line in f:
		if ("@<TRIPOS>BOND") in line:
			break
		arr = line.split()
		charge += float(arr[8])
print('{:.2f}'.format(charge))
