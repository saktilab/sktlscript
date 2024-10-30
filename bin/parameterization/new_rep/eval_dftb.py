#!/usr/bin/env python3

import subprocess
import tempfile
import os
import scipy.constants
import numpy

def Bool2YesNo(value):
    if (value):
        return "Yes"
    else:
        return "No"

class Skfileinfo:
    def __init__(self, json: map, basepath = None):
        self.json = json      
        if basepath is not None:
            abspath = os.path.abspath(basepath)
            skfiles_node = self.json.get('slaterkosterfiles', None)
            if (skfiles_node):
                if ('type2filenames' in skfiles_node):
                    inputfile = skfiles_node['type2filenames']['prefix']
                    if (not os.path.isabs(inputfile)):
                        skfiles_node['type2filenames']['prefix'] = os.path.join(abspath, inputfile)
                else:
                    for k, v in self.json.items():
                        if (not os.path.isabs(v)):
                            skfiles_node[k] = os.path.join(abspath, inputfile)
        
    def write(self, fout):
        skfiles_node = self.json.get('slaterkosterfiles', None)
        print('\n+Hamiltonian = +DFTB{', file=fout)
        if ('type2filenames' in skfiles_node):
            node = skfiles_node['type2filenames']
            print('    !SlaterKosterFiles = !Type2Filenames{', file=fout)
            print('      !Prefix = "{}"'.format(node['prefix']), file=fout)
            print('      !Suffix = "{}"'.format(node['suffix']), file=fout)
            print('      !Separator = "{}"'.format(node['separator']), file=fout)
            print('      !LowerCaseTypeName = {}'.format(Bool2YesNo(node['lowercasetypename'])), file=fout)
            print('    }', file=fout)

            print('    /MaxAngularMomentum = {', file=fout)
            for k, v in self.json['maxangularmomentum'].items():
                print('      {} = "{}"'.format(k, v), file=fout)
            print('    }', file=fout)

        else:
            print('    !SlaterKosterFiles = {', file=fout)
            for k, v in skfiles_node.items():
                print('      {} = {}'.format(k, v), file=fout)
            print('    }', file=fout)

            print('    /MaxAngularMomentum = {', file=fout)
            for k, v in self.json['maxangularmomentum'].items():
                print('      {} = "{}"'.format(k, v), file=fout)
            print('    }', file=fout)
        print('}', file=fout)

def evaluate_atom(dftb_path, dftb_input, skfileinfo):
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, 'dftb_in.hsd'), 'w') as fout:
            #copy the template input
            with open(dftb_input, 'r') as fin:
                fout.write(fin.read())
            #add SK file info
            skfileinfo.write(fout)
        # with open(os.path.join(tmp, 'dftb_in.hsd'), 'r') as f:
        #     print(f.read())
        res = subprocess.run([dftb_path], stdout=subprocess.PIPE, cwd=tmp)
        if (res.returncode == 0):
            with open(os.path.join(tmp,'detailed.out'), 'r') as fin:
                for line in fin:
                    if 'Total energy' in line:
                        return float(line.split()[2])*scipy.constants.physical_constants['atomic unit of electric potential'][0]
        else:
            raise RuntimeError(res.stdout.decode())            

def write_system(system, fout):
    geometry = system['geometry']
    #write geometry
    symbols = []
    coords = []
    for sym, coord in geometry['coordinates']:
        symbols.append(sym)
        coords.append(coord)
    elements = list(set(symbols))
    pbc = False
    if ('lattice_vectors' in geometry):
        pbc = True
    print('!Geometry = !GenFormat{', file=fout)
    # fractional = False
    if (not pbc):
        print('{} C'.format(len(symbols)), file=fout)
    else:
        if (geometry['fractional_coordinates']):
            print('{} F'.format(len(symbols)), file=fout)
            # fractional = True
        else:
            print('{} S'.format(len(symbols)), file=fout)
    print(' '.join(elements),file=fout)
    for i in range(len(symbols)):
        print('{} {} {:20.12F} {:20.12F} {:20.12F}'.format(
            i+1, elements.index(symbols[i])+1, *coords[i]
        ), file=fout)
    if (pbc):
        print("0.0 0.0 0.0", file=fout)
        for i in range(3):
            print('{:20.12F} {:20.12F} {:20.12F}'.format(
             *geometry['lattice_vectors'][i]
        ), file=fout)
    print('}', file=fout)

    #write Hamiltonian
    print('\n+Hamiltonian = +DFTB{', file=fout)
    print('   !Charge = {}'.format(geometry['charge']), file=fout)
    if (pbc):
        if ('supercell_folding' in geometry['kpoints']):
            print('   !KPointsAndWeights = SupercellFolding{', file=fout)
            for line in geometry['kpoints']['supercell_folding']:
                print('    {} {} {}'.format(*line), file=fout)
            print('    }', file=fout)
        else:
            raise RuntimeError('Not implemented')
    print('}', file=fout)



def evaluate_system(dftb_path, system, skfileinfo):
    force_factor = 51.42208619083232
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, 'dftb_in.hsd'), 'w') as fout:
            #copy the template input
            with open(system['template_input'], 'r') as fin:
                fout.write(fin.read())
            #add SK file info
            skfileinfo.write(fout)
            write_system(system, fout)
            print('*Options = {', file=fout)
            print('    !WriteAutotestTag = Yes', file=fout)
            print('}', file=fout)
            print('*Analysis = {', file=fout)
            print('    !CalculateForces = Yes', file=fout)
            print('}', file=fout)
            
        # with open(os.path.join(tmp, 'dftb_in.hsd'), 'r') as f:
            # print(f.read())
        res = subprocess.run([dftb_path], stdout=subprocess.PIPE, cwd=tmp)
        if (res.returncode == 0):
            # print(res.stdout.decode())
            with open(os.path.join(tmp,'autotest.tag'), 'r') as fin:
                for line in fin:
                    if line.startswith('total_energy'):
                        line = next(fin)
                        total_energy = float(line)*scipy.constants.physical_constants['atomic unit of electric potential'][0]
                    elif line.startswith('repulsive_energy'):
                        line = next(fin)
                        rep_energy = float(line)*scipy.constants.physical_constants['atomic unit of electric potential'][0]
                    elif line.startswith('repulsive_forces'):
                        entry = line.split(':')[3]
                        nrow = int(entry.split(',')[1])
                        rep_forces = []
                        for i in range(nrow):
                            line = next(fin)
                            rep_forces.append(list(map(float, line.split())))
                    elif line.startswith('forces'):
                        entry = line.split(':')[3]
                        nrow = int(entry.split(',')[1])
                        forces = []
                        for i in range(nrow):
                            line = next(fin)
                            forces.append(list(map(float, line.split())))


                elec_energy = (total_energy - rep_energy)
                elec_forces = ((numpy.array(forces)-numpy.array(rep_forces))*force_factor).tolist()
                return total_energy, elec_energy, elec_forces
        else:
            raise RuntimeError(res.stdout.decode())  


def evaluate(eval_settings, ref_systems):

    skfileinfo = Skfileinfo(eval_settings['skfileinfo'], os.curdir)

    for atom_symbol, atom_energy in eval_settings['atomic_energy'].items():
        if (type(atom_energy) != float):
            res = evaluate_atom(eval_settings['dftb_path'], atom_energy, skfileinfo)
            eval_settings['atomic_energy'][atom_symbol] = res

    results = {}
    for ind, system in enumerate(ref_systems['system']):
        print('Evaluating {}/{}'.format(ind+1, len(ref_systems['system'])))
        _, elec_energy, elec_forces = evaluate_system(eval_settings['dftb_path'],
            system, skfileinfo
        )
        # print(energy, forces)
        results[system['uuid']] = {}
        results[system['uuid']]['elec_energy'] = elec_energy
        results[system['uuid']]['elec_forces'] = elec_forces

        cohesive_energy = 0.0
        
        geometry = system['geometry']    
        for sym, _ in geometry['coordinates']:
            cohesive_energy += eval_settings['atomic_energy'][sym]

        cohesive_energy -= elec_energy
        cohesive_energy -= system['cohesive_energy']

        target_rep_force = numpy.array(system['forces']) - numpy.array(elec_forces)


        results[system['uuid']]['target_rep_energy'] = cohesive_energy
        results[system['uuid']]['target_rep_forces'] = target_rep_force.tolist()


    return results