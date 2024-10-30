#!/usr/bin/env python3

import pymatgen as mg 
import pymatgen.symmetry.analyzer as sa
import sys

struct = mg.Structure.from_file(sys.argv[1])

sga = sa.SpacegroupAnalyzer(struct)
# pg = sa.PointGroupAnalyzer(struct, tolerance=0.3, eigen_tolerance=0.01, matrix_tol=0.1)
print(sga.get_space_group_operations())

# print(sga.get_symmetrized_structure())
symmops = sga.get_symmetry_operations()
so = sa.SpacegroupOperations("Pmmn",59,symmops)


# print(so.are_symmetrically_equivalent(sites1=list(struct.sites[0]),sites2=list(struct.sites[1]),symm_prec=0.001))
# for site in struct.sites:
#     print(site)
    # print(so.are_symmetrically_equivalent(site,site, symm_prec=0.001))
