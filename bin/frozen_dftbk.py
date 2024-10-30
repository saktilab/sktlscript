#!/usr/bin/env python3
import numpy as np 
import sys

f = open(sys.argv[1], "r")

f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
moveatoms = []
index = 0
for line in f:
    arr = line.split()
    index += 1
    if (arr[3] == "F" and arr[4] == "F" and arr[5] == "F"):
        moveatoms.append(index)

for i in moveatoms:
    print(i, end=' ')        
#if(sys.argv[2]=='commas'):
#    for i in moveatoms:
#        print(i, end=',')
#else:
#    for i in moveatoms:
#        print(i, end=' ')
