#!/usr/bin/env python3

import sys
import math
import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description='Create NVE input from NVT MD')

parser.add_argument('output', metavar='output', type=str)
parser.add_argument('traject', metavar='traject', type=str)
parser.add_argument('velocity', metavar='velocity', type=str)
parser.add_argument('-l', '--latt', type=str, required=True, help='Lattice')

parser.add_argument('-t', '--temperature', type=float, default=300.0 )
parser.add_argument('-i', '--template_input', type=str, required=True )
parser.add_argument('-s', '--skip', type=int, default=0 )

opt = parser.parse_args(sys.argv[1:])


def get_coord_at_step(filename, step):
    with open(opt.traject, 'r') as f:
        for line in f:
            nat = int(line)
            title = next(f)    
            arr = title.split()
            cur_step = int(arr[9])
            if (cur_step == step):
                coords = []
                for i in range(nat):
                    line = next(f)
                    coords.append(line.strip())
                return coords
            else:
                for i in range(nat):
                    next(f)

def get_veloc_at_step(filename, step):
    with open(opt.velocity, 'r') as f:
        for line in f:
            nat = int(line)
            title = next(f)    
            arr = title.split()
            cur_step = int(arr[9])
            if (cur_step == step):
                coords = []
                for i in range(nat):
                    line = next(f)
                    coords.append(line.split()[1:4])
                return coords
            else:
                for i in range(nat):
                    next(f)

print_interval = 1
with open(opt.latt, 'r') as f:
    for line in f:
        if 'MD' in line and 'PRINT' in line:
            m = re.search('PRINT=([0-9]*)', line)
            print_interval = int(m.groups()[0])
            

TV = []
SK_section = []
with open(opt.latt, 'r') as f:
    sec_ind = 1
    for line in f:
        if (sec_ind == 3):
            SK_section.append(line.rstrip())
        if (len(line.strip()) == 0):
            sec_ind += 1
        if ('TV' in line):
            TV.append(line.strip())
good_step = -1
good_temp = 0
nat = 0
spin = 1
charge = 0

good_info = []
with open(opt.output, 'r') as f:
    for line in f:
        if 'STEP NO.' in line:
            arr = line.split()
            step = int(arr[9])


            next(f)
            title = next(f)
            temp = float(title.split()[2])

            if (step % print_interval != 0 ) :
                continue
            if (temp > opt.temperature*0.99 and temp < opt.temperature*1.01):
                good_info.append((step, temp)) 
#                good_step = step
#                good_temp = temp
        elif 'Total number of atoms' in line :
            nat = int(line.split()[5])
        elif 'Spin multiplicity' in line:
            spin = int(line.split()[3])
        elif 'Charge of system' in line:
            charge = int(line.split()[4])

nlines = len(good_info)
good_step, good_temp = good_info[nlines-1-opt.skip]
print(good_step, good_temp)

coords = get_coord_at_step(opt.traject, good_step)
velocity = get_veloc_at_step(opt.velocity, good_step)


with open('veloc.dat', 'w') as fout:
    for line in velocity:
        print (' '.join(line), file=fout)

section_index = 1
with open(opt.template_input, 'r') as f:
    with open('dftb.inp', 'w') as fout:
        for line in f:
            if (line.strip() == ''):
                section_index += 1
                if (section_index == 3):
                    print('',file=fout)
                    break
            print(line,end='', file=fout)
        print("\n".join(SK_section), file=fout)
        print('{} {} {}'.format(nat, charge, spin), file=fout)
        for line in coords:
            print(line, file=fout)
        for line in TV: 
            print(line, file=fout)
        print('', file=fout)
