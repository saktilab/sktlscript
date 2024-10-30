#!/usr/bin/env python3

import sys
from itertools import islice

with open(sys.argv[1],'r') as f:
    nAtom = int(next(f))

with open(sys.argv[1],'r') as f:
    for line in f:
        if line.rstrip() == "energy":
            print(list(islice(f, nAtom))[-1])
            break
