#!/usr/bin/env python3

import numpy as np
import sys


step = 1
with open(sys.argv[1], "r") as f:
    for line in f:
        nAtom = int(line)
        comments = next(f)
        arr = comments.split()
        step += 1

        


