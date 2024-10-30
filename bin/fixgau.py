#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as f:
	symb = []
	x = []
	y = []
	z = []
	nAtom = next(f)
	next(f)
	for line in f:
		arr = line.split()
		symb.append(arr[0])
		x.append(arr[1])
		y.append(arr[2])
		z.append(arr[3])

frozen_atoms = sys.argv[2:]
froze_label = []
for i, ele in enumerate(symb):
	if str(i+1) in frozen_atoms:
		froze_label.append('-1')
	else:
		froze_label.append('0')

with open("geom_fix.xyz",'w') as fout:
	for ele, label, x,y,z in zip(symb,froze_label,x,y,z):
		print("{} {} {} {} {}".format(ele,label,x,y,z),file=fout)