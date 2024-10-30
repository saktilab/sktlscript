#!/usr/bin/env python3
import sys

symb = []
x =[]
y =[]
z =[]

with open(sys.argv[1],'r') as f:
    next(f)
    next(f)
    for line in f:
        arr = line.split()
        symb.append(arr[0])
        x.append(arr[1])
        y.append(arr[2])
        z.append(arr[3])

filename=sys.argv[1].split('.')[0]

NoAtom = {'C':'6','Ni':'28','P':'15','N':'7','H':'1','O':'8','S':'16','Cl':'17'}

with open('{}.gms'.format(filename),'w') as fout:
    for symb, x, y, z in zip(symb,x,y,z):
        print("{} {} {} {} {}".format(symb, float(NoAtom[symb]),x,y,z),file=fout)
