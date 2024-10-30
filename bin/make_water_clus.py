#!/usr/bin/env python3

import sys
import pymatgen as mg
import scipy.constants
import subprocess
import io
import os
import shutil


target_nat = sys.argv[2]

solvent_file = sys.argv[1]
solvent_name_array = os.path.splitext(os.path.basename(solvent_file))
solvent_path = os.path.abspath(os.path.dirname(solvent_file))
solvent_name = solvent_name_array[0]
solvent_ext = solvent_name_array[-1]

script_path = os.path.dirname(os.path.abspath(__file__))


if (solvent_ext != ".xyz"):
    raise ValueError('I need a .xyz file for input')

print('Solvent name: {}'.format(solvent_name))

target_density = 1.0
#target_density = float(sys.argv[2])
if solvent_name in density_map:
    target_density = density_map[solvent_name]
target_ion_con = 1.0

Avo = scipy.constants.physical_constants['Avogadro constant'][0]
packmol = os.path.expanduser("~/bin/packmol")

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

radius = (mass*nmole*1.0E24/(target_density*Avo))**(1.0/3.0)
n_ions = int(target_ion_con*box_size**3.0*1.0E-27*Avo+0.5)

print('Target # of atoms in system: {}'.format(target_nat))
print('Target system density: {}'.format(target_density))
print('Solvent size: {}'.format(nat))
print('# of solvent molecules to match target # of atoms: {}'.format(nmole))
print('Size of the simulation box: {:20.12F}'.format(box_size))
print('# of Lithium cations to match the target ion concentration of {:8.4} M: {}'.format(target_ion_con,n_ions))


if not os.path.exists('prep'):

    print('Making pure solvent system using packmol')
    os.makedirs('prep', exist_ok=True)
    os.chdir('prep')

    shutil.copyfile('{}.pdb'.format(os.path.join(solvent_path, solvent_name)), '{}.pdb'.format(solvent_name))

    with open('packmol.inp', 'w') as fout:
        print("""
        tolerance 2.0
        filetype pdb
        output system_init.pdb 

        structure {}.pdb 
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


    shutil.copyfile(os.path.join(script_path, 'make_box_pure.in'), "make_box.in")
    shutil.copyfile('{}.lib'.format(os.path.join(solvent_path, solvent_name)), 'sol.lib')
    shutil.copyfile('{}.frcmod'.format(os.path.join(solvent_path, solvent_name)), 'sol.frcmod')
    result = run_command("tleap", ["-f", "make_box.in"])
    print(result)


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

