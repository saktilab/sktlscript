#!/usr/bin/env python3
import sys
import statistics as st
import numpy as np
f = open(sys.argv[1], "r")
fout = open("nonzeros.dat", "w")
data = []
for line in f:
    if(float(line) != 0.0):
    	data.append(float(line))

for i in data:
	print(i, file=fout)
average = st.mean(data)
median = st.median(data)
minimum = min(data)
maximum = max(data)
data = np.array(data, dtype="float")
std = np.std(data)
print("{} from {} data".format(average,len(data)))
print("Minimum: {}".format(minimum))
print("Maximum: {}".format(maximum))
print("Std: {}".format(std))
print("Median: {}".format(median))
