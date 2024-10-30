#!/usr/bin/env python3

import numpy as np
from sklearn.linear_model import LinearRegression
import argparse
import sys

parser = argparse.ArgumentParser('Calculate Regression from arbritary Data in x, y format')
parser.add_argument('-i','--input',type=str, help='Data file')

opt = parser.parse_args(sys.argv[1:])

x = []
y = []

with open(opt.input, 'r') as f:
    next(f)
    for line in f:
        arr = line.split()
        x.append(float(arr[0]))
        y.append(float(arr[1]))

x = np.array(x).reshape((-1,1))
y = np.array(y)

model = LinearRegression()
model.fit(x,y)

r_sq = model.score(x,y)

print('Rsquare value = {}'.format(r_sq))

print('slope = {}'.format(model.coef_))

print('intercept = {}'.format(model.intercept_))

