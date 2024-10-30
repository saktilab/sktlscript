#!/usr/bin/env python3
import sys
import pymatgen as mg
import random
import math
import pymatgen.transformations.standard_transformations as trans
import argparse

N_H2O = 0
H2O_DISTANCE = 3.0
HOF_DISTANCE = 3.0


parser = argparse.ArgumentParser(description='Pack H2O in HOF')
parser.add_argument('-n', default=0, type=int)
parser.add_argument('filename', metavar='filename', type=str)
opt = parser.parse_args(sys.argv[1:])

if (opt.n == 0):
    parser.print_usage()
    sys.exit(1)


N_H2O = opt.n

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


hof_str = mg.Structure.from_file(opt.filename)
h2o_sites = []

current_num_h2o = 0
while current_num_h2o < N_H2O:
    pt = [random.uniform(0.0,1.0) for i in range(0, 3)]
    o_coord = hof_str.lattice.get_cartesian_coords(pt)
    o_site = mg.PeriodicSite("O", o_coord, hof_str.lattice,coords_are_cartesian=True).to_unit_cell

    # check if it's too close to other water
    valid = True
    for site in h2o_sites:
        if (o_site.distance(site) < H2O_DISTANCE):
            # print ("Collision")
            valid = False
            break

    #check if the carbon is too close to other atoms in mof_str
    neighbors = hof_str.get_neighbors(o_site, HOF_DISTANCE)
    if len(neighbors) > 0:
        valid = False

    if (valid):
        current_num_h2o += 1
        print ("Add H2O #{}".format(current_num_h2o))
        o_coord = o_site.coords.copy()
        h_coord_1 = o_site.coords.copy()
        h_coord_1[2] += 1.3
        h_coord_2 = o_site.coords.copy()
        h_coord_2[2] -= 1.3

        temp_mole = mg.Molecule(["O", "H", "H"], [o_coord, h_coord_1, h_coord_2])
        com = temp_mole.center_of_mass
        com_mole = temp_mole.get_centered_molecule()
        h2o_mole = rotate(com_mole)
        h2o_mole.translate_sites(range(0,len(h2o_mole.sites)), com)

        for site in h2o_mole.sites:
            h2o_site = mg.PeriodicSite(site.specie, site.coords, lattice = hof_str.lattice, coords_are_cartesian=True).to_unit_cell
            h2o_sites.append(h2o_site)

h2o_sites.sort(key=lambda x : x.specie)
for site in h2o_sites:
    new_site = site.to_unit_cell
    hof_str.sites.append(new_site)

# mof_str.sort()
hof_str.to(fmt='poscar', filename='H2O_{}.POSCAR'.format(N_H2O))