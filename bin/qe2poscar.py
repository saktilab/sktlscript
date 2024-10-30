#!/usr/bin/env python
from ase.io import read,write
import os
import sys

#os.system("sed -i 's/magn=//' {}".format(sys.argv[1]))
struct = read(sys.argv[1],format='espresso-out')

fname = sys.argv[1].split('.')[0]

write('{}.vasp'.format(fname),struct)
