#!/usr/bin/env python3
import sys

fout=open("serial_O.dat", "w")
fout2=open("serial_H.dat", "w")
step=0
with open(sys.argv[1], "r") as f:
    natom = int(f.readline())
    next(f)
    step +=1
    iatom = 0
    init_O = []
    init_H = []
    for line in f:
        iatom +=1
        if (iatom < natom):
            arr = line.split()
            symbol = arr[0]
            if (symbol == "O"):
                init_O.append(iatom)
            if (symbol == "H"):
                init_H.append(iatom)
    for i in init_O:
        print(i,end=' ', file=fout)

    for i in init_H:
        print(i, end=' ',file=fout2)
