#!/usr/bin/env python3
import sys
import argparse
import os
from shutil import copyfile

parser = argparse.ArgumentParser(description='### A protocol to automatically run NVE, NVT, and NPT simulations in NAMD ###')
parser.add_argument('-pdb', '--pdb', type=str, help='NAMD pdb file')
parser.add_argument('-psf', '--psf', type=str, help='NAMD psf file')
parser.add_argument('-par', '--par', type=str, help='NAMD parameter file')
parser.add_argument('-tcl', '--tcl', type=str, help='NAMD config file')
opt = parser.parse_args(sys.argv[1:])

###First of first is NVE simulations
def run_nve():
    protocoler = ['minimization','heating', 'equilibration', 'production']
    input_file = [opt.pdb, opt.par, opt.psf, opt.tcl]
    if not os.path.exists('NVE'):
        os.mkdir('NVE')
        for protocol in protocoler:
            os.mkdir(protocol)
    for file in input_file:
        copyfile(file, "./NVE/minimization/")
    os.chdir('NVE/minimization/')
run_nve()

