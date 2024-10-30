#!/usr/bin/env python3
import sys
import numpy as np

data = []

with open(sys.argv[1],"r") as f:
    for line in f:
        arr = line.split()
        data.append(arr[0:])

data = np.array(data, dtype='float')
np.savetxt("table_{}".format(sys.argv[1]), data,fmt='%.2f', delimiter=' &', newline='\\\\\n')
