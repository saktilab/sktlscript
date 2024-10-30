#!/usr/bin/env python3
import sys

with open(sys.argv[1], 'r') as f:
    next(f)
    atom = []
    for line in f:
        if "ATOM" in line:
            arr = line.split()
            atom.append(arr[5])
serial = 1
print(serial)
for x,y in zip(atom[::],atom[1::]):
    if x == y:
        serial = serial
    else:
        serial += 1
    print(serial)