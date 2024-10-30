#!/usr/bin/env python3
import pymatgen as mg
import argparse 
import sys
import os
import pymatgen.io.vasp as pv
import pymatgen.analysis.defects.core  as dc
from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.defects.generators import VacancyGenerator,Vacancy

#Remove annoying warnings!
import warnings

parser = argparse.ArgumentParser(description='###Build vacancy structures###')

warnings.filterwarnings("ignore")
parser.add_argument('POSCAR', metavar='POSCAR', type=str)
parser.add_argument('-a', '--super_a', type=int, default=1, help='length of supercell in a axis')
parser.add_argument('-b', '--super_b', type=int, default=1, help='length of supercell in b axis')
parser.add_argument('-c', '--super_c', type=int, default=1, help='length of supercell in c axis')

opt = parser.parse_args(sys.argv[1:])
bulk = mg.Structure.from_file(opt.POSCAR)
ala = VacancyGenerator(bulk,include_bv_charge=True)
ala.__getattribute__
for i in range(len(ala.structure.sites)):
    vac = Vacancy(bulk,defect_site=ala.structure.sites[i])
    vac = vac.generate_defect_structure()
    vac.to(filename="vac_{}.vasp".format(i), fmt="poscar")
# ala = Vacancy(bulk,charge=-2,defect_site=ala.structure.sites[32])
# ala = ala.generate_defect_structure()
# ala.to(filename="ala.vasp", fmt="poscar")



    