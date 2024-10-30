#!/usr/bin/env python3

import pymatgen as mg
import pymatgen.core.units

def get_energy_in_au(item):
    if (item.unit.lower() == 'h'):
        return item.energy
    elif (item.unit.lower() == 'ev'):
        return mg.core.units.Energy(item.energy, 'eV').to('Ha')
    elif (item.unit.lower() == 'kcal/mol'):
        return item.energy/627.50947415
    else:
        raise ValueError('Unknown unit: {}'.format(item.unit))
    
def get_length_in_au(item):
    if (item.unit.lower() == 'bohr'):
        return item.distance
    elif (item.unit.lower() == 'aa'):
        return item.distance*mg.core.units.ang_to_bohr
    else:
        raise ValueError('Unknown unit')
    
def get_energy_from_au(value, item):
    if (item.unit.lower() == 'h'):
        return value
    elif (item.unit.lower() == 'ev'):
        return mg.core.units.Energy(value, 'Ha').to('eV')
    elif (item.unit.lower() == 'kcal/mol'):
        return value*627.50947415
    else:
        raise ValueError('Unknown unit: {}'.format(item.unit))