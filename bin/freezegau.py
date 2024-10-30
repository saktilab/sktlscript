#!/usr/bin/env python3
### This is a script to convert the serial of frozen atoms into the Gaussian coordinate format. The script automatically assign '-1' string to the frozen atoms.
import sys

symb = []
x =[]
y =[]
z =[]

with open(sys.argv[1],'r') as f:
    natom = int(next(f))
    next(f)
    for line in f:
        arr = line.split()
        symb.append(arr[0])
        x.append(arr[1])
        y.append(arr[2])
        z.append(arr[3])

filename=sys.argv[1].split('.')[0]

freezeatoms = sys.argv[2].split(',')

freezeindex = []
for i in range(natom):
    if str(i+1) in freezeatoms:
        freezeindex.append('-1')
    else:
        freezeindex.append('0')

with open('{}.gau'.format(filename),'w') as fout:
    for symb, freezeindex, x, y, z in zip(symb,freezeindex,x,y,z):
        print("{} {} {} {} {}".format(symb,freezeindex,x,y,z),file=fout)
