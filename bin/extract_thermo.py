#!/usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(description="A simple script to extract thermochemical information from Gaussian 09 output")
parser.add_argument("-i","--input", type=str, default="inp.log",help="Gaussian output file name")

opt = parser.parse_args(sys.argv[1:])

with open(opt.input, 'r') as f:
    for line in f:
        if "Sum of electronic and zero-point Energies" in line:
            arr = line.split()
            print("Electronic energy (E_elec) = {} hartree".format(float(arr[6])))
        if "Sum of electronic and thermal Energies=" in line:
            arr = line.split()
            print("Internal energy (E) = {} hartree".format(float(arr[6])))
        if "Sum of electronic and thermal Enthalpies=" in line:
            arr = line.split()
            print("Enthalpy (H) = {} hartree".format(float(arr[6])))
        if "Sum of electronic and thermal Free Energies=" in line:
            arr = line.split()
            print("Gibbs free energy (G) = {} hartree".format(float(arr[7])))
