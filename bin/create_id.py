#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as f:
    next(f)
    id = []
    for line in f:
        arr = line.split()
        id.append(int(arr[2]))
    for n, i in enumerate(id):
        if i==-1:
            id[n] = id[n-1]
        #if i==152:
         #   id[n] = 153
with open("ala.dat", 'w') as fout:
    for i in id:
        print("index {}".format(i), file=fout)
