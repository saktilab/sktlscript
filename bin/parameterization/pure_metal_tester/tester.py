#!/usr/bin/env python3

import sys
import os
import argparse
import yaml
import pprint
import subprocess
import tempfile
import shutil
import glob
import pymatgen
import evtesting
import distutils.dir_util
import pymatgen.io.vasp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
import eos
import math
from utils import *
import elastic
import surface
import scipy.constants


dftbplus_path = 'dftb+'
template_dftb_input_file = os.path.abspath('inputs/dftb_in.hsd')
template_surface_input_file = os.path.abspath('inputs/surface_dftb2.hsd')
template_atom_input_file = os.path.abspath('inputs/atom.hsd')
skfile_path = os.path.abspath('skfiles/')

DFT_reference = {
    'optimized_files': {
        'fcc': 'testing/fcc.vasp',
        'bcc': 'testing/bcc.vasp',
        'sc': 'testing/sc.vasp',
        'hcp': 'testing/hcp.vasp',
    },
    'optimized_geoms': {},
}

ev_testing = {
    'points': 11,
    'strain': 0.05    
}

Hartree2Ev = scipy.constants.physical_constants['Hartree energy in eV'][0]
eV2Joule = scipy.constants.physical_constants['electron volt-joule relationship'][0]

def print_testing_info():
    for k, v in DFT_reference['optimized_files'].items():
        print ('Loading reference geometry for {} phase: {}'.format(k, v))
        with open(v, 'r') as f:
            struct = pymatgen.Structure.from_str(f.read(), fmt='poscar')
        DFT_reference['optimized_geoms'][k] = sga(struct).get_conventional_standard_structure()
    
print_testing_info()

       
class TestingFailure(BaseException):
    def __init__(self, reason):
        self.reason = reason


def test_ev_curve(root_path, phase, dirn, tolerence=0.02):
    if (phase in DFT_reference['optimized_geoms']):
        print('='*80)
        print('Testing E-V curves of {} phase...'.format(phase))
        optimized_struct = DFT_reference['optimized_geoms'][phase]
        deformed_structs = evtesting.make_strain_structures(optimized_struct, dirn=dirn, mdr=ev_testing['strain'], NoP=ev_testing['points'])
        print('------------------------------------------------')
        print('Testing {} E-V'.format(dirn))
        print('------------------------------------------------')
        print('{:>6s} {:>20s} {:>20s}'.format('Strain', 'Volume(A^3)/Atom', 'Energy(H)/Atom'))
        print('------------------------------------------------')
        strain_energies = []
        energies = []
        volumes = []
        for strain, struct_ in deformed_structs:
            running_path = os.path.join(root_path, 'ev', phase, dirn, '{:.3f}'.format(strain))
            res = evaluate_dftb(struct_, path=running_path, skfile_path=skfile_path, temp_input_path=template_dftb_input_file, dftbplus_path=dftbplus_path, geom_opt=False, opt_lattice=False)
            energy = res['energy_final']
            energy_per_atom = energy/len(struct_.sites)
            volume = struct_.volume
            volume_per_atom = volume/len(struct_.sites)
            print ('{:6.3f} {:20.12f} {:20.12f}'.format(strain, volume_per_atom, energy_per_atom))
            volumes.append(volume/(0.52977**3.0))
            strain_energies.append((strain, energy_per_atom))
            energies.append(energy)
        print('------------------------------------------------')

        if (dirn == 'VOL'):
            # for v, e in zip(volumes, energies):
            #     print(v, e)
            eos_res = eos.calculate_eos(volumes, energies)
            print('BM EOS fitting:')
            print('V0: {:10.4f}, B0: {:10.4f}, BP: {:10.4f}, log(chi): {:8.3f}'.format(
                eos_res['V0']*0.529177**3.0/len(struct_.sites), eos_res['B0'], eos_res['BP'], eos_res['chi']
                )
            )

            if (math.isnan(eos_res['chi'])):
                print('BM EOS fitting Failed.')
                raise TestingFailure('EV {}'.format(dirn))


        strain_energies.sort(key=lambda x: x[1])
        if (abs(strain_energies[0][0]) > tolerence ):
            print('The minimal of EV curve is located at {:.2f}. Failed.'.format(strain_energies[0][0]))
            raise TestingFailure('EV {}'.format(dirn))
        else:
            print('The minimal of EV curve is located at {:.2f}. Passed.'.format(strain_energies[0][0]))


def run_ev_testing(root_path, force=True):
    print('Running testing for E-V curves...')
    tests = [
        ('fcc','VOL', 0.02),
        ('bcc','VOL', 0.02),
        ( 'sc','VOL', 0.03),
        ('hcp','VOL', 0.02),
        ('hcp','COA', 0.02),
        ('hcp','BOA', 0.02),
    ]
    for phase, dirn, tol in tests:
        try:
            test_ev_curve(root_path, phase, dirn, tol)
            sys.stdout.flush()
        except TestingFailure  as ex:
            if (force==False):
                pass
            else:
                raise ex
    print(' All passed')
    print('='*80)



def run_geometry_optimization_testing(root_path, force=True):
    print('='*80)
    print('Testing geometry optimization for all phases...')
    phase_energies = []

    res = {}
    res['geom'] = {}
    
    for phase, struct in DFT_reference['optimized_geoms'].items():
        try:
            running_path = os.path.join(root_path, 'opt', phase)
            
            sgb4 = sga(struct).get_space_group_number()

            print('Phase: {}'.format(phase))
            struct_11 = struct.copy()
            struct_11.scale_lattice(struct_11.volume*1.1)
            res_11 = evaluate_dftb(struct_11, running_path, skfile_path=skfile_path, temp_input_path=template_dftb_input_file, dftbplus_path=dftbplus_path, geom_opt=True, opt_lattice=True, opt_fixangles=True)
            sgaf_11 = sga(res_11['geom_final']).get_space_group_number()

            struct_09 = struct.copy()
            struct_09.scale_lattice(struct_09.volume*0.9)
            res_09 = evaluate_dftb(struct_09, running_path, skfile_path=skfile_path, temp_input_path=template_dftb_input_file, dftbplus_path=dftbplus_path, geom_opt=True, opt_lattice=True, opt_fixangles=True)
            sgaf_09 = sga(res_09['geom_final']).get_space_group_number()

            print('Space group: before {}, after(11) {}, after(09) {}'.format(sgb4, sgaf_11, sgaf_09 ))
            
            print('          {:>10s} {:>10s} {:>10s} {:>7s} {:>7s} {:>7s} {:>10s}'.format(
                'a','b','c','alpha','beta','gamma','volume'))
            print('-'*80)
            b4array = list(struct.lattice.abc) + list(struct.lattice.angles) + [struct.volume]
            print('Before    : {:10.4f} {:10.4f} {:10.4f} {:7.2f} {:7.2f} {:7.2f} {:10.4f}'.format(
                *b4array
            ))
            afarray_11 = list(res_11['geom_final'].lattice.abc) + list(res_11['geom_final'].lattice.angles) + [res_11['geom_final'].volume]
            print('After(11) : {:10.4f} {:10.4f} {:10.4f} {:7.2f} {:7.2f} {:7.2f} {:10.4f}'.format(
                *afarray_11
            ))
            afarray_09 = list(res_09['geom_final'].lattice.abc) + list(res_09['geom_final'].lattice.angles) + [res_09['geom_final'].volume]
            print('After(09) : {:10.4f} {:10.4f} {:10.4f} {:7.2f} {:7.2f} {:7.2f} {:10.4f}'.format(
                *afarray_09
            ))
            diff_11 = [(y-x)/x*100 for x, y in zip(b4array,afarray_11)]
            print('Delta11[%]: {:10.2f} {:10.2f} {:10.2f} {:7.2f} {:7.2f} {:7.2f} {:10.2f}'.format(
                *diff_11
            ))
            diff_09 = [(y-x)/x*100 for x, y in zip(b4array,afarray_09)]
            print('Delta09[%]: {:10.2f} {:10.2f} {:10.2f} {:7.2f} {:7.2f} {:7.2f} {:10.2f}'.format(
                *diff_09
            ))
            print('-'*80)

            if (sgb4 != sgaf_11 or sgb4 != sgaf_09):
                raise TestingFailure('Geom. opt. failed for {} phase'.format(phase))

            max_lc_error = 5.0
            max_volume_error = 10.0
            if (max([abs(x) for x in diff_11[:-1]]) > max_lc_error or abs(diff_11[-1])>max_volume_error or 
                max([abs(x) for x in diff_09[:-1]]) > max_lc_error or abs(diff_09[-1])>max_volume_error ) :
                raise TestingFailure('Geom. opt. failed for {} phase'.format(phase))

            phase_energies.append((phase, min(res_09['energy_final'], res_11['energy_final'])/float(len(struct.sites))))
            
            if (res_09['energy_final']<res_11['energy_final']) :
                res['geom'][phase] = res_09['geom_final']
            else:
                res['geom'][phase] = res_11['geom_final']
        except TestingFailure  as ex:
            print(ex)
            if (force==False):
                pass
            else:
                raise ex
        sys.stdout.flush()

    phase_energies.sort(key=lambda x: x[1])
    min_ene = phase_energies[0][1]
    
    print('-'*80)
    running_path = os.path.join(root_path, 'atomic')
    atom_res = evaluate_dftb_input(running_path, skfile_path=skfile_path, temp_input_path=template_atom_input_file, dftbplus_path=dftbplus_path)
    atom_energy = atom_res['energy_final'] 
    print('Atomic energy: {} H, {} eV'.format(atom_energy, atom_energy*Hartree2Ev))
    print('Cohesive energies:')
    ind = 0
    for phase, energy in phase_energies:
        if (ind == 0):
            print('{:5s} {:+10.4f} eV'.format(phase, (energy-atom_energy)*Hartree2Ev))
        else:
            print('{:5s} {:+10.4f} eV'.format(phase, (energy-min_ene)*Hartree2Ev))
        ind += 1


    print('='*80)
    res['phase_diff_energy'] = phase_energies
    return res

def run_elastic_constant_testing(root_path, phase, optimized_struct):
    print('Computing elestic constants for phase {}'.format(phase))
    evaluator = elastic.ElaticEvaluator()
    res = evaluator.compute_elastic_constants(root_path, optimized_struct, skfile_path=skfile_path, temp_input_path=template_dftb_input_file, dftbplus_path=dftbplus_path)
    # print(res)

def calc_surface_energy(root_path, phase, optimized_struct, temp_dftb_input=None):
    evaluator = surface.SurfaceEnergyEvaluator()
    res = evaluator.compute_surface_energy(root_path, optimized_struct, skfile_path=skfile_path, temp_input_path=template_surface_input_file, dftbplus_path=dftbplus_path)

curdir = os.path.abspath(os.curdir)
with DummyBlock('./temp') as tmpdir: #tempfile.TemporaryDirectory() as tmpdir_:

    print('Temp dir: {}'.format(tmpdir))
    os.chdir(tmpdir)
        
    
#    run_ev_testing(tmpdir, force=False)
    
    # Geom Opt. Testing
    # Check optimized Space Group
    optimized_structures = run_geometry_optimization_testing(tmpdir, force=False)
    
    lowest_phase = optimized_structures['phase_diff_energy'][0][0]
    # Elastic Constant Test
    # for phase, struct in optimized_structures['geom'].items():
#    run_elastic_constant_testing(tmpdir, lowest_phase, optimized_structures['geom'][lowest_phase])

    # Surface energy
    # calc_surface_energy(tmpdir, lowest_phase, optimized_structures['geom'][lowest_phase])

    
os.chdir(curdir)






