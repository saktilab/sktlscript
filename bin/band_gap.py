#!/usr/bin/env python3
import sys



with open(sys.argv[1], "r") as f:
    for line in f:
        arr = line.split()
        if float(arr[1]) == 0:
            print(arr[0])

with open(sys.argv[2], "r") as f:
    for line in f:
        arr = line.split()

