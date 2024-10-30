#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import scipy
import peakutils 


f = open(sys.argv[1], "r")
x = []
y = []

next(f)
next(f)
for line in f:
    arr = line.split()
    x.append(arr[0])
    y.append(float(arr[1]))


indices = peakutils.indexes(y, thres=0.6, min_dist=0.1)
x_peaks = []
y_peaks = []
for i in indices:
    x_peaks.append(x[i])
    y_peaks.append(y[i])

x_peaks_ind = np.argsort(y_peaks)
x_peaks_ind[:] = x_peaks_ind[::-1]
for i in x_peaks_ind:
    print(x_peaks[i])

