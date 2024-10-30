#!/usr/bin/env python
import sys
import random
import math
from pymatgen.core import Structure, Molecule, PeriodicSite
import pymatgen.transformations.standard_transformations as trans
import argparse

N_GAS = 0



parser = argparse.ArgumentParser(description='Pack molecule in framework')
parser.add_argument('-n', default=0, type=int)
parser.add_argument('-m', '--mole', type=str)
parser.add_argument('-dg', default=2.0, type=float,help='Gas intermolecular distance')
parser.add_argument('-df', default=2.0, type=float,help='Gas to frame intermolecular distance')
parser.add_argument('-f', '--frame', type=str)
opt = parser.parse_args(sys.argv[1:])

GAS_DISTANCE = opt.dg
FRAME_DISTANCE = opt.df

if (opt.n == 0):
    parser.print_usage()
    sys.exit(1)


N_GAS = opt.n

def distance(pt1, pt2):
    return math.sqrt(sum((pt1-pt2)**2))

def rotate(mole):
    angles = [random.uniform(0.0,360.0) for i in range(0, 3)]
    op1 = trans.RotationTransformation([1,0,0], angles[0])    
    op2 = trans.RotationTransformation([0,1,0], angles[1])
    op3 = trans.RotationTransformation([0,0,1], angles[2])
    old_mole = mole
    new_mole = op1.apply_transformation(old_mole)
    old_mole = new_mole
    new_mole = op2.apply_transformation(old_mole)
    old_mole = new_mole
    new_mole = op3.apply_transformation(old_mole)
    return new_mole

frame_str = Structure.from_file(opt.frame)
mole_str = Molecule.from_file(opt.mole)
mole = mole_str.get_centered_molecule()
gas_sites = []

current_num_gas = 0
while current_num_gas < N_GAS:
    pt = [random.uniform(0.0,1.0) for i in range(0, 3)]
    cm_coord = frame_str.lattice.get_cartesian_coords(pt)
    new_mole = mole.copy()
    new_mole.translate_sites(vector=cm_coord)
    new_mole = rotate(new_mole)

    valid = True
    temp_sites = []
    for site in new_mole:
        new_site = PeriodicSite(site.specie, site.coords, frame_str.lattice, coords_are_cartesian=True).to_unit_cell      
    
        # check if it's too close to other GAS    
        for osite in gas_sites:
            if (new_site().distance(osite) < GAS_DISTANCE):
                # print ("Collision")
                valid = False
                break
        if (valid == False):
            break          
        
        #check if the carbon is too close to other atoms in frame_str
        neighbors = frame_str.get_neighbors(new_site(), FRAME_DISTANCE)
        if len(neighbors) > 0:
            valid = False
        
        temp_sites.append(new_site())
    
    if (valid):
        current_num_gas += 1 
        print ("Add Gas #{}".format(current_num_gas))
        gas_sites += temp_sites

gas_sites.sort(key=lambda x : x.specie)           
for site in gas_sites:
    new_site = site.to_unit_cell
    frame_str.sites.append(new_site())

# frame_str.sort()
frame_str.to(fmt='POSCAR', filename='GAS_{}.vasp'.format(N_GAS))

