#!/usr/bin/env python3
import pymatgen as mg
import pymatgen.analysis.adsorption as pa
import pymatgen.core.structure as st
import pymatgen.core.surface as sf
import argparse 
import numpy as np
import sys
import warnings
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser(description='Build possible adsorbate-adsorbent interactions for the surface reaction')
parser.add_argument('-s', '--slab', type=str, help='structure POSCAR of the surface' )
parser.add_argument('-height','--height', type=float, help='height criteria for selection of surface sites' )
parser.add_argument('-d','--distance', type=float, help='distance from the coordinating ensemble of atoms along the miller index for the site (i. e. the distance from the slab itself)' )
parser.add_argument('-pi','--put_inside', type=bool, help='whether to put the site inside the cell ' )
parser.add_argument('-sr','--symm_reduce', type=float, help='symmetry reduction threshold' )
parser.add_argument('-nr','--near_reduce', type=float, help='near reduction threshold' )
parser.add_argument('-noh','--no_obtuse_hollow', type=bool, help='flag to indicate whether to include obtuse triangular ensembles in hollow sites' )
# parser.add_argument('-l','--list', type=bool, default=True, help='List the number of possible sites' )
# parser.add_argument('-p','--print_to_file', type=bool, default=False, help='print the structure to VASP input')
parser.add_argument('-sel','--selective_dynamics', type=bool, default=False, help='create the selective dynamics options for VASP calculations')
parser.add_argument('-t','--type', type=str, help='Type of adsorption site based on the list' )
parser.add_argument('-in','--index', type=int, help='Site index' )
parser.add_argument('-ad','--ad', type=str, help='Structure of adsorbate from file' )


opt = parser.parse_args(sys.argv[1:])

########ANIMATION#########
import itertools
import threading
import time
import sys

done = False
#here is the animation
def animate():
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done:
            break
        sys.stdout.write('\rloading ' + c)
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write('\rDone!     ')

t = threading.Thread(target=animate)
t.start()

#long process here
time.sleep(10)
done = True

slab_coord = mg.Structure.from_file(opt.slab)
adsorbent = pa.AdsorbateSiteFinder(slab_coord, selective_dynamics=opt.selective_dynamics, height=opt.height)
sites = adsorbent.find_adsorption_sites(distance=opt.distance, put_inside=opt.put_inside, symm_reduce=opt.symm_reduce, near_reduce=opt.near_reduce, positions=['ontop', 'bridge', 'hollow'], no_obtuse_hollow=opt.no_obtuse_hollow)

#Get all sites
hollow = sites.get('hollow')
bridge = sites.get('bridge')
ontop = sites.get('ontop')

#Defining the adsorbate structure from file
ad = st.Molecule.from_file(opt.ad)

print('###############THE LIST OF ADSORPTION SITES####################')
print('hollow:{}, bridge:{}, ontop:{}'.format(len(hollow), len(bridge), len(ontop)))

if (opt.type == 'hollow'):    
    adsorbed = adsorbent.add_adsorbate(ad, hollow[opt.index], repeat=None, reorient=True)
    adsorbed.to(filename='adsorbed_{}_{}.vasp'.format(opt.type, opt.index), fmt='poscar')    

if (opt.type == "bridge"):
    adsorbed = adsorbent.add_adsorbate(ad, bridge[opt.index], repeat=None, reorient=True)
    adsorbed.to(filename='adsorbed_{}_{}.vasp'.format(opt.type, opt.index), fmt='poscar')

if (opt.type == "ontop"):
    adsorbed = adsorbent.add_adsorbate(ad, ontop[opt.index], repeat=None, reorient=True)
    adsorbed.to(filename='adsorbed_{}_{}.vasp'.format(opt.type, opt.index), fmt='poscar')



