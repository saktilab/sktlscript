#!/usr/bin/env python3

import sys
import pymatgen as mg
import scipy.constants
import subprocess
import io
import math
import os
import shutil
from psfgen import PsfGen


# density_map = {
# "APN_curved": 0.951,
# "APN_linear": 0.951,
# "ATN"       : 0.786,
# "BC"        : 1.14,
# "DEC"       : 0.975,
# "DIOX"      : 1.034,
# "DME"       : 0.867,
# "DMF"       : 0.944,
# "DMI"       : 1.056,
# "DMSO"      : 1.10,
# "EC"        : 1.321,
# "EMC"       : 1.006,
# "FEC"       : 1.41,
# "GBL"       : 1.12,
# "GN_curved" : 0.995 ,
# "GN_linear" : 0.995, 
# "GVL"       : 1.05,
# "MAN"       : 0.956 ,
# "MPN"       : 0.937,
# "MSL"       : 1.200 ,
# "M-THF"     : 0.86,
# "NE"        : 1.045,
# "NM"        : 1.127,
# "NMO"       : 1.2   ,
# "NMP"       : 1.028,
# "PC"        : 1.204,
# "SFC"       : 1.261,
# "THF"       : 0.889,
# "TMP"       : 1.197 ,
# "VC"        : 1.355
# }

solvent_file = sys.argv[1]
solvent_name_array = os.path.splitext(os.path.basename(solvent_file))
solvent_path = os.path.abspath(os.path.dirname(solvent_file))
solvent_name = solvent_name_array[0]
solvent_ext = solvent_name_array[-1]

script_path = os.path.dirname(os.path.abspath(__file__))


if (solvent_ext != ".xyz"):
    raise ValueError('I need a .xyz file for input')

print('Solvent name: {}'.format(solvent_name))

target_density = 1
#target_density = float(sys.argv[2])
# if solvent_name in density_map:
#     target_density = density_map[solvent_name]
# target_ion_con = 1.0
target_nat = int(sys.argv[2])

Avo = scipy.constants.physical_constants['Avogadro constant'][0]
packmol = os.path.expanduser("~/bin/packmol_bin")

def run_command(command, args=[], input_str=None):
    command_list = [command]
    if (len(args) > 0):
        command_list.extend(args)

    if (input_str is not None):
        res =  subprocess.run(command_list, stdout=subprocess.PIPE, input=input_str.encode(encoding='ascii'))
    else:
        res =  subprocess.run(command_list, stdout=subprocess.PIPE)
    return res.stdout.decode()


mole = mg.Molecule.from_file(solvent_file)
mass = sum([x.specie.atomic_mass for x in mole.sites])
nat = len(mole.sites)

nmole = int(target_nat / nat)

box_size = (mass*nmole*1.0E24/(target_density*Avo))**(1.0/3.0)
# radius = (3*mass*nmole*1.0e24/(4*math.pi*target_density*Avo))**(1.0/3.0)
# n_ions = int(target_ion_con*box_size**3.0*1.0E-27*Avo+0.5)

print('Target # of atoms in system: {}'.format(target_nat))
print('Target system density: {}'.format(target_density))
print('Number of atom in one water molecule: {}'.format(nat))
print('# of water molecules to match target # of atoms: {}'.format(nmole))
print('Size of the simulation water sphere: {:20.12F}'.format(box_size))
# print('# of Lithium cations to match the target ion concentration of {:8.4} M: {}'.format(target_ion_con,n_ions))


if not os.path.exists('prep'):

    print('Making pure solvent system using packmol')
    os.makedirs('prep', exist_ok=True)
    os.chdir('prep')

    with open('packmol.inp', 'w') as fout:
        print("""
        tolerance 2.0
        filetype xyz
        output system_init.xyz 

        structure {}.xyz
          number {}
          inside cube -{} -{} -{} {}
        end structure 

        """.format(solvent_name, nmole, box_size/2.0, box_size/2.0, box_size/2.0, box_size)
        , file=fout)

    print('='*80)
    with open('packmol.inp', 'r') as f:
        print(f.read())
    print('='*80)

    result = run_command(packmol)
    with open('packmol.out', 'w') as f:
        print(result, file=f)
    os.system('~/bin/xyz2tip4p_2005.py system_init.xyz > system_init.pdb')
    gen = PsfGen() 
    gen.read_topology('{}.rtf'.format(os.path.join(solvent_path, solvent_name)))
    gen.add_segment(segid='P0', pdbfile='system_init.pdb') 
    gen.read_coords(segid='P0', filename='system_init.pdb')
    gen.regenerate_angles() 
    gen.regenerate_dihedrals() 
    gen.write_psf(filename='system.psf')
    gen.write_pdb(filename='system.pdb')

    shutil.copyfile('{}.par'.format(os.path.join(solvent_path, solvent_name)), 'SOL.par')

    os.chdir('../')

os.makedirs('nve/mini', exist_ok=True)
os.makedirs('nve/heat', exist_ok=True)
os.makedirs('nve/equil', exist_ok=True)
os.makedirs('nve/prod', exist_ok=True)
os.makedirs('nvt', exist_ok=True)
os.makedirs('npt', exist_ok=True)

with open(os.path.join(script_path, 'minimization.conf'), 'r') as fin, open(os.path.join('nve','mini', 'namd_in.conf'), 'w') as fout:
    for line in fin:
        print(line, file=fout, end='')
    print('cellBasisVector1 {} 0 0'.format(box_size), file=fout)
    print('cellBasisVector2 0 {} 0'.format(box_size), file=fout)
    print('cellBasisVector3 0 0 {}'.format(box_size), file=fout)

shutil.copyfile(os.path.join(script_path, 'heating.conf'), os.path.join('nve', 'heat', 'namd_in.conf'))
shutil.copyfile(os.path.join(script_path, 'equil.conf'), os.path.join('nve', 'equil', 'namd_in.conf'))
shutil.copyfile(os.path.join(script_path, 'nve.conf'), os.path.join('nve', 'prod', 'namd_in.conf'))
shutil.copyfile(os.path.join(script_path, 'nvt.conf'), os.path.join('nvt', 'namd_in.conf'))
shutil.copyfile(os.path.join(script_path, 'npt.conf'), os.path.join('npt', 'namd_in.conf'))
shutil.copyfile(os.path.join(script_path, 'npt_rest.conf'), os.path.join('npt', 'namd_rest.conf'))

# os.makedirs('npt', exist_ok=True)
# radpot = radius + 2.0
# with open(os.path.join(script_path, 'minimization.conf'), 'r') as fin, open(os.path.join('nve','mini', 'namd_in.conf'), 'w') as fout:
#     for line in fin:
#         print(line, file=fout, end='')
#     print('sphericalBCr1 {}'.format(radpot), file=fout)

# with open(os.path.join(script_path, 'heating.conf'), 'r') as fin, open(os.path.join('nve','heat', 'namd_in.conf'), 'w') as fout:
#     for line in fin:
#         print(line, file=fout, end='')
#     print('sphericalBCr1 {}'.format(radpot), file=fout)

# with open(os.path.join(script_path, 'equil.conf'), 'r') as fin, open(os.path.join('nve','equil', 'namd_in.conf'), 'w') as fout:
#     for line in fin:
#         print(line, file=fout, end='')
#     print('sphericalBCr1 {}'.format(radpot), file=fout)

# with open(os.path.join(script_path, 'nve.conf'), 'r') as fin, open(os.path.join('nve','prod', 'namd_in.conf'), 'w') as fout:
#     for line in fin:
#         print(line, file=fout, end='')
#     print('sphericalBCr1 {}'.format(radpot), file=fout)

# with open(os.path.join(script_path, 'nvt.conf'), 'r') as fin, open(os.path.join('nvt', 'namd_in.conf'), 'w') as fout:
#     for line in fin:
#         print(line, file=fout, end='')
#     print('sphericalBCr1 {}'.format(radpot), file=fout)
#     print('set temperature {}'.format(float(sys.argv[3])), file=fout)
#     print('temperature $temperature', file=fout)
#     print('langevinTemp $temperature', file=fout)
    
# shutil.copyfile(os.path.join(script_path, 'heating.conf'), os.path.join('nve', 'heat', 'namd_in.conf'))
# shutil.copyfile(os.path.join(script_path, 'equil.conf'), os.path.join('nve', 'equil', 'namd_in.conf'))
# shutil.copyfile(os.path.join(script_path, 'nve.conf'), os.path.join('nve', 'prod', 'namd_in.conf'))
# shutil.copyfile(os.path.join(script_path, 'nvt.conf'), os.path.join('nvt', 'namd_in.conf'))
# shutil.copyfile(os.path.join(script_path, 'npt.conf'), os.path.join('npt', 'namd_in.conf'))
# shutil.copyfile(os.path.join(script_path, 'npt_rest.conf'), os.path.join('npt', 'namd_rest.conf'))

