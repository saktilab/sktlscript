#!/usr/bin/env python3

import sys


NAtom = 0

with open(sys.argv[1], 'r') as f:
	for line in f:
		if 'MOLECULE' in line:
			next(f)
			NAtom += int(next(f).split()[0])
		if 'ATOM' in line:
			charge = []
			for i in range(NAtom):
				charge.append(float(next(f).split()[8])) 

print(sum(charge))
