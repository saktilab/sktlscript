#!/usr/bin/env python3

import pymatgen as mg
import numpy.linalg as la
import numpy as np
import math
import sys


ori_str = mg.Structure.from_file(sys.argv[1])
ori_str.make_supercell([2,2,1])


cart_coords = ori_str.cart_coords

v_a, v_b, v_c = ori_str.lattice.matrix

norm_b_prime = np.cross(v_a, v_c)
norm_b_prime /= la.norm(norm_b_prime)

if (np.sign(norm_b_prime[0]) != np.sign(v_b[0])):
    norm_b_prime = -norm_b_prime

new_v_b = la.norm(v_b)*math.cos(30.0/180.0*math.pi)*norm_b_prime



new_lv = [v_a, new_v_b, v_c]

new_str = mg.Structure(new_lv, ori_str.species, cart_coords, coords_are_cartesian = True)
new_str.to(fmt='poscar', filename='POSCAR.rhom.super')

final_cart_coords = []
final_species = []
for i in range(len(new_str.sites)):
        new_site = new_str.sites[i].to_unit_cell
        if (new_site.frac_coords[0] < 0.5):
            final_cart_coords.append(new_site.coords)
            final_species.append(new_site.specie)

new_v_a = v_a/2.0
final_lv = [ new_v_a, new_v_b, v_c ]

final_str = mg.Structure(final_lv, final_species, final_cart_coords, coords_are_cartesian = True)
final_str.to(fmt='poscar', filename='POSCAR.rhom')
