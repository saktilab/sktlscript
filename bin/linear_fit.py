#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import statsmodels.api as sm
import matplotlib
matplotlib.use('PS')
matplotlib.rc('text', usetex=True)
#matplotlib.use('Qt5agg')
import matplotlib.pyplot as plt
import warnings

#Remove annoying warnings!
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description='General linear fitting scheme')
parser.add_argument('data', metavar='data', type=str)
parser.add_argument('-s', '--start', type=int, default=0, help='StartStep')
opt = parser.parse_args(sys.argv[1:])


data_file = opt.data

x = []
ys = []
ncol = -1
lineno = 0
headers = []
with open(data_file, 'r') as f:
    for line in f:
        lineno +=1
        if (line.startswith('#')):
            if (lineno == 1):
                headers = line.split()[1:]
            continue
        arr = line.split()
        x.append(float(arr[0]))
        ys.append(list(map(float, arr[1:])))
        
        if (ncol == -1):
            ncol = len(arr) -1
        else:
            if (ncol != (len(arr) -1)):
                raise ValueError('The number of columns of data file is not consistent!')
        

print('There are {} columns of data.'.format(ncol))
print()

print(headers)
if (len(headers) != ncol):
    headers = []
    for col in range(ncol):
        headers.append('Column {}'.format(col+1))

for col in range(ncol):
    print ('Column {}:'.format(col+1))
    
    all_data_x = np.array(x, dtype=np.double).reshape(-1, 1)
    data_x = all_data_x[opt.start:]
    
    all_data_y = [d[col] for d in ys]
    data_y = np.array(all_data_y[opt.start:], dtype=np.double).reshape(-1, 1)
    X = sm.add_constant(data_x)
    ols = sm.OLS(data_y, X)
    ols_result = ols.fit()
    
    summ = ols_result.summary()
    print(summ)
    
    print()
    
    predicted_y = ols_result.predict(X)

    plt.plot(x, all_data_y, label=headers[col])
    plt.plot(data_x, predicted_y, 'k--')
    
# plt.legend(loc=2)
# plt.xlabel('Time (ps)')    
# plt.ylabel('MSD (Angstrom$^2$)')
# plt.savefig('msd.pdf', dpi=1200, format='pdf')
# plt.show()    



