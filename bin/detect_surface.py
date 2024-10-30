#!/usr/bin/env python3
import sys

# symbol = []
# x = []
# y = []
# z = []
# index = []
index = -1
surf_ele=[]
with open(sys.argv[1], "r") as f:
    next(f)
    next(f)
    next(f)
    next(f)
    next(f)
    next(f)
    next(f)
    next(f)
    for line in f:
        arr=line.split()
        symbol = arr[3]
        z = float(arr[2])
        index += 1
        if (symbol == str(sys.argv[2]) and z > float(sys.argv[3])):
            surf_ele.append(index)
            with open("surface_elements.dat","w") as fout:
                for index in surf_ele:
                    print(index, end=' ', file=fout)
print("Number of surface elements: {}".format(len(surf_ele)))
