#!/usr/bin/env python3
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
opt = parser.parse_args(sys.argv[1:])
struct = mg.Structure.from_file(opt.POSCAR)
