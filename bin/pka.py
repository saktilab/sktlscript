#!/usr/bin/env python3

import sys

dG = float(sys.argv[1])
T = float(sys.argv[2])

R = 8.314

pKa = dG/(2.303*R*T)*1000
print(pKa)