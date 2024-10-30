#!/usr/bin/env python3
import sys

f1 = open(sys.argv[1],"r")
f2 = open(sys.argv[2], "r")

time = []
h1 = []
h2 = []
next(f1)
for line in f1:
    arr = line.split()
    time.append(float(arr[0]))
    h1.append(int(arr[1]))
next(f2)
for line in f2:
    arr = line.split()
    h2.append(int(arr[1]))
htot = [x + y for x,y in zip(h1, h2)]

with open("htot.dat", "w") as fout:
    for i,j in zip(time,htot):
        print(i, j, file=fout)