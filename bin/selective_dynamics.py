#!/usr/bin/env python3
import numpy as np 
import sys
import os

f = open(sys.argv[1], "r")

f.readline()
f.readline()
tv1 = f.readline().split()
tv2 = f.readline().split()
tv3 = f.readline().split()
f.readline()
f.readline()
f.readline()
f.readline()
moveatoms = []
frozen = []
index = 0
for line in f:
    arr = line.split()
    index += 1
    if (arr[3] == "T" and arr[4] == "T" and arr[5] == "T"):
        moveatoms.append(index)
    else:
        frozen.append(index)

with open("move_atoms.dat", 'w') as fout:
    for i in moveatoms:
        print(i, end=' ', file=fout)

with open("frozen_atoms.dat", 'w') as fout:
    for i in frozen:
        print(i, end=' ', file=fout)

os.system("~/bin/frozen.sh {}".format(sys.argv[2]))

with open("dftb.inp", 'w') as fout:
    print("""SCC=(MAXITER=200 )
OPT=(MAXITER=1000 GCONV=1e-3)
DISP=(DISPTYPE=5)
DC=FALSE

Optimasi geometri

    3
    O 2
     O-O.skf O-Zn.skf O-C.skf
    Zn 3
     Zn-O.skf Zn-Zn.skf Zn-C.skf
    C 2
     C-O.skf C-Zn.skf C-C.skf
        """, file=fout)
    natoms = len(moveatoms) + len(frozen)
    print(natoms, "0", "1", end='\n', file=fout)
    with open("geom.xyz", 'r') as f:
        for line in f:
            print(line, end='', file=fout)
    print("TV", tv1[0], tv1[1], tv1[2], file=fout)
    print("TV", tv2[0], tv2[1], tv2[2], file=fout)
    print("TV", tv3[0], tv3[1], tv3[2], file=fout)
    print('', file=fout)