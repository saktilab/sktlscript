#!/usr/bin/env python3

from ase.io import read, write

import argparse 
import sys
import re

parser = argparse.ArgumentParser(description='Convert different format of molecular structure')
parser.add_argument('-i', '--input', type=str, help='Input geometry' )
parser.add_argument('-if','--input_format', type=str, help='Input geometry format' )
parser.add_argument('-o', '--output', type=str, help='Output geometry' )
parser.add_argument('-of','--output_format', type=str, help='Output geometry format')
parser.add_argument('-s', '--sort', type=bool, help='Sort the element')
parser.add_argument('-e', '--element', type=str, help='Element list to be sorted')


opt = parser.parse_args(sys.argv[1:])

input_structure = read(opt.input, format=opt.input_format)
output_structure = write(opt.output, input_structure, format=opt.output_format)



if (opt.sort == True):
    element = re.findall('[A-Z][^A-Z]*', opt.element)
    input_structure = input_structure.sort(key=lambda x: element.index(str(x.specie))) 


if(opt.output_format == 'vasp'):
    f = open(opt.output, "r")
    names = f.readline().split() 
    f.close()
    f = open(opt.output, "r")
    contents = f.readlines()
    f.close()
    contents.insert(5, ' '.join(names))
    contents.insert(6, ' '.join("\n"))
    f = open(opt.output, "w")        
    f.writelines(str(x) for x in contents)
    f.close()




