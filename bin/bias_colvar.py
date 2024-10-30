#!/usr/bin/env python3
import sys

gau_height = []
colvar = []
steps = []
step = 0
with open(sys.argv[1], "r") as f:
    for line in f:
        if "Gaussian height" in line:
            arr = line.split()
            gau_height.append(arr[3])
        if "Coordinate" in line:
            arr = line.split()
            colvar.append(arr[2])
        step += 1
        steps.append(step)
for steps, x, y in zip(steps, colvar,gau_height):
    print(steps,x,y)