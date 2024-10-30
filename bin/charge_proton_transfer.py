#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import statsmodels.api as sm

parser = argparse.ArgumentParser(description='Compute Grotthuss Proton transfer')
parser.add_argument('time_charge', metavar='time_charge', type=str)
opt = parser.parse_args(sys.argv[1:])

times = []
max_charge_atoms = []
step = 0
with open(opt.time_charge, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        arr = line.split()
        times.append(float(arr[0]))
        if (float(arr[2]) < -0.65):
            max_charge_atoms.append(max_charge_atoms[-1])
        else:
            max_charge_atoms.append(int(arr[1]))


num_pros = 0

iList = max_charge_atoms[0]
jList = max_charge_atoms[0]
if (max_charge_atoms[1]!= max_charge_atoms[0]):
    num_pros = 1
    iList = max_charge_atoms[1]

lList = max_charge_atoms[0]
time = times[0]

data_x = []
data_y = []

data_x.append(time)
data_y.append(num_pros)
print ('#{:>9s} {:>8s} {:>8s} '.format('Time(ps)', '#Proton', 'MaxIndex'))
print ('{:10.4f} {:>8d} {:>8d} '.format(time, num_pros, max_charge_atoms[1]))
for i in range(2, len(max_charge_atoms)):
    time = times[i]
    if (max_charge_atoms[i] == max_charge_atoms[i-1]):
        pass
    else:
        if (max_charge_atoms[i] != iList and max_charge_atoms[i] == jList):
            if (iList == lList):
                num_pros += 1
            else:
                num_pros -= 1
                lList = max_charge_atoms[i]
        else:
            num_pros += 1
        jList = iList
        iList = max_charge_atoms[i]

    data_x.append(time)
    data_y.append(num_pros)
    print ('{:10.4f} {:>8d} {:>8d}'.format(time, num_pros, max_charge_atoms[i]))

    
# all_data_x = np.array(data_x, dtype=np.double).reshape(-1, 1)

# all_data_y = np.array(data_y, dtype=np.double).reshape(-1, 1)

# ols = sm.OLS(all_data_y, all_data_x)
# ols_result = ols.fit()

# summ = ols_result.summary()
# print(summ)
