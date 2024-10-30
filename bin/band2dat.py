#!/usr/bin/env python3

import sys

bandout = sys.argv[1]
shift_energy = float(sys.argv[2])

bands_alpha = []
bands_beta  = []

with open(bandout, 'r') as f:
    for line in f:
        if ('KPT' in line):
            bands = bands_alpha
            if ('SPIN' in line):
                if (int(line.split()[3])==2):
                    bands = bands_beta
            band = []
            while True:
                line = next(f)
                if (len(line.strip())==0):
                    break
                band.append(float(line.split()[0])-shift_energy)
            bands.append(band)

with open('bands_alpha.dat', 'w') as f:
    bands = bands_alpha
    for i, line in enumerate(bands):
        formatstr = '{:12.6f}'*len(bands[0])
        final_line = formatstr.format(*line)
        print("{} {}".format(i, final_line), file=f)

with open('bands_beta.dat', 'w') as f:
    bands = bands_beta
    for i, line in enumerate(bands):
        formatstr = '{:12.6f}'*len(bands[0])
        final_line = formatstr.format(*line)
        print("{} {}".format(i, final_line), file=f)
