#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import sys
import os

#x = []
y = []
data = open(sys.argv[1], 'r')
for line in data:
    arr = line.split()
 #   x.append(arr[0])
    y.append(arr[0])
#x = np.array(x, dtype='float')
x = np.arange(int(len(y)))
y = np.array(y, dtype='float')
yhat = savgol_filter(y,int(sys.argv[2]),3)
c = np.column_stack([x,yhat])
filename, file_extension = os.path.splitext(sys.argv[1])
np.savetxt("smooth_{}.dat".format(filename), c)
    
# plt.plot(x,y)
# plt.plot(x, yhat, color='red')
# plt.show()
