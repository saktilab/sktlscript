#!/usr/bin/env python3

from scipy.constants import k, h
from math import log
import sys

R = 8.314
r = float(sys.argv[1])
T = float(sys.argv[2])
Ea = R*T*log(k*T/h/r) # in J/mol
joule2kcal = 0.000239006
Ea = Ea*joule2kcal

print("PT barrier = {} kcal/mol".format(Ea))