#!/usr/bin/env python3

import sys
f = open(sys.argv[1], "r")
indices = []
index = -1
for line in f:
    index += 1
    if ('*' in line):
        indices.append(index)
indices = [int(x) for x in indices]
print(indices)