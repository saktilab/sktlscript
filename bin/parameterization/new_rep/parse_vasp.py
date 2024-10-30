#!/usr/bin/env python3

import sys
import pymatgen.io.vasp
import xml.etree.ElementTree as ET
import json
import numpy
import uuid
import slugify
import os

def reversed_string(a_string):
    return a_string[::-1]

def parse_system_charge(vasprun, doc):
    NELECT = vasprun.parameters['NELECT']
    elem_atominfo = doc.find('atominfo')
    nums = 0
    for a in elem_atominfo.findall("array"):
        if a.attrib["name"] == "atomtypes":
            nums = sum([float(rc.findall("c")[0].text)*float(rc.findall("c")[3].text) for rc in a.find("set")])
    return int(nums-NELECT)
    
def parse_force_stress(doc):
    elem_calc = doc.find('calculation')
    forces = []
    stress = []
    for elem in elem_calc.findall('varray'):
        if (elem.attrib['name'] == 'forces'):
            forces = [[float(i) for i in v.text.split()] for v in elem]
        if (elem.attrib['name'] == 'stress'):
            stress = [[float(i) for i in v.text.split()] for v in elem]
    return forces, stress

def StructureToJson(struct, kpoints, charge):
    res = {}
    coords = []
    for site in struct.sites:
        coords.append( [str(site.specie), list(site.frac_coords)] )
    res['fractional_coordinates'] = True
    res['coordinates'] = coords
    node_kpoints = {}

    kpts = numpy.zeros((3,3), dtype=numpy.int).tolist()
    for i in range(3):
        kpts[i][i] = kpoints.kpts[0][i]
        
    
    if (kpoints.style ==  pymatgen.io.vasp.Kpoints.supported_modes.Gamma):
        kpts.append([0.0, 0.0, 0.0])
    elif (kpoints.style ==  pymatgen.io.vasp.Kpoints.supported_modes.Monkhorst):
        kpts.append([0.5, 0.5, 0.5])

    node_kpoints['supercell_folding'] = kpts
    res['kpoints'] = node_kpoints
    res['lattice_vectors'] = struct.lattice.matrix.tolist()
    res['scaling_factor'] = 1.0
    res['charge'] = charge
    return res

def compute_cohesive_energy(total_energy, struct, atomic_energies):
    all_needed_atoms = [str(x.specie) for x in struct.sites]
    energy = 0.0
    for atom in all_needed_atoms:
        if atom not in atomic_energies:
            raise ValueError('Not enough atomic energy')
    for site in struct.sites:
        energy += atomic_energies[str(site.specie)]
    return energy - total_energy

class VaspResult:
    def __init__(self, vasprun_file, atomic_energies={}):
        self.vasprun = pymatgen.io.vasp.Vasprun(vasprun_file, parse_potcar_file=False, parse_eigen=False, parse_projected_eigen=False)

        doc = ET.parse(vasprun_file)
        self.forces, self.stress = parse_force_stress(doc)
        self.charge = parse_system_charge(self.vasprun, doc)
        with open(vasprun_file, 'r') as f:
            self.uuid = str(uuid.uuid5(uuid.NAMESPACE_X500, f.read()))
        self.name = reversed_string(slugify.slugify(reversed_string(os.path.dirname(os.path.abspath(vasprun_file))), max_length=20, word_boundary=True))
        self.atomic_energies = atomic_energies
    def to_dict(self):
        res = {}
        res['method'] = self.vasprun.parameters.get('GGA', 'unknown')
        res['basis'] = '/'.join([ x['titel'] for x in self.vasprun.potcar_spec])
        res['geometry'] = StructureToJson(self.vasprun.final_structure, self.vasprun.kpoints, self.charge)
        res['forces'] = self.forces
        res['total_energy'] = self.vasprun.final_energy
        try:
            res['cohesive_energy'] = compute_cohesive_energy(self.vasprun.final_energy, self.vasprun.final_structure, self.atomic_energies)
        except:
            pass

        res['name'] = self.name
        res['template_input'] = 'dftbinp/crystal_energy.hsd'
        res['uuid'] = self.uuid
        return res
