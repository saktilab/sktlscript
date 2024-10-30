#!/usr/bin/env python3

import sys
import pymatgen as mg
import pandas as pd
import statistics
import math
import functools 
import itertools
import json
from collections import namedtuple


class PythonObjectEncoder(json.JSONEncoder):
    def default(self, obj):
        print(type(obj))
        return json.JSONEncoder.default(self, obj)
        

def get_poscar_frame(filename):
    with open(filename, 'r') as f:
        while 1:
            lines = []
            line = f.readline()
            if not line.strip():
                StopIteration()
                break;
            lines.append(line.strip())

            for _ in range(5):
                line = next(f)
                lines.append(line.strip())
            
            line = next(f)
            lines.append(line.strip())
            nat = sum(map(int,line.split()))
            line = next(f)
            lines.append(line.strip())
            for _ in range(nat):
                line = next(f)
                lines.append(line.strip())
            yield '\n'.join(lines)

class MOFStructureInfo:
    def find_metals(self):
        self.metal_indexes = [x for x in range(len(self.struct.sites)) 
        if self.struct.sites[x].specie.is_alkali or 
           self.struct.sites[x].specie.is_alkaline or 
           self.struct.sites[x].specie.is_transition_metal or 
           self.struct.sites[x].specie.is_rare_earth_metal
           ]
    
    
    def __init__(self, struct):
        self.struct = struct
        self.find_metals()
        self.get_bond_length_angle_indexes()
        
    


    def as_dict(self):
        d = {
            'geometry': self.struct.as_dict(),
            'coordination_number': self.coordination_number,
            'bond_length_indexes': self.bond_length_indexes,
            'bond_angle_indexes': self.bond_angle_indexes
        }
        return d
    
    def get_bond_length_angle_indexes(self, max_dis = 2.4):
        pair_index = []
        angle_index = []
        cn = []
        for i in self.metal_indexes:
            metal = self.struct.sites[i]
            neighbors = self.struct.get_neighbors(metal, max_dis, include_index=True)
            to_remove = [i for i, x in enumerate(neighbors) if str(x[0].specie)=='H' and x[1] > 1.4]
            for index in reversed(to_remove): 
                del neighbors[index]
            cn.append(len(neighbors))
            neighbor_indexes = [x[2] for x in neighbors]
            pair_index += [[i, int(x[2])] for x in neighbors]
            angle_index += [[int(x[0]), i, int(x[1])] for x in itertools.combinations(neighbor_indexes, 2)]


        pair_index += [[m1, m2] for m1, m2 in itertools.combinations(self.metal_indexes, 2) if self.struct.get_distance(m1, m2) < 4.0]
            

        self.coordination_number = cn
        self.bond_length_indexes = pair_index
        self.bond_angle_indexes = angle_index
        return 

class ComparisonErrorModel:
    def __init__(self, name):
        self.name = name
        self.data = {}    
       
    def get_stat_value(self):
        res = {}
        
        all_data = []
        detailed = {}
        for key, tuples in self.data.items():
            data = [x[-1] for x in tuples]
            
            points = len(data)
            res_item = {}
            res_item['npoints'] = points
            res_item['max'] = max(map(math.fabs, data))
            res_item['mad'] = statistics.mean(map(math.fabs, data))
            res_item['mse'] = statistics.mean(data)
            res_item['rms'] = math.sqrt(statistics.mean([x**2 for x in data]))
                      
            all_data += data
            detailed[key] = res_item

        res['npoints'] = len(all_data)
        res['max'] = max(map(math.fabs,all_data))
        res['mad'] = statistics.mean(map(math.fabs, all_data))
        res['mse'] = statistics.mean(all_data)
        res['rms'] = math.sqrt(statistics.mean([x**2 for x in all_data]))
        res["detailed"] = detailed
        return res

    def as_dict(self):
        d = {
            'stats' : self.get_stat_value(),
            'data': self.data
        }
        return d


class BondLengthComparisonErrorModel(ComparisonErrorModel):
    def __init__(self, name):
        super().__init__(name)
        self.DataType = namedtuple('DataEntry', ['atom1', 'a1', 'atom2', 'a2', 'ref', 'calc', 'error'])
        
    
    def add_value(self, atom1, index1, atom2, index2, ref_value, calc_value):
        pair = '{}-{}'.format(*list(sorted([atom1, atom2])))
        error = calc_value - ref_value
        if (pair in self.data):
            self.data[pair] += [[atom1, index1, atom2, index2,
            ref_value, calc_value, error]]
        else:
            self.data[pair] = [[atom1, index1, atom2, index2, 
            ref_value, calc_value, error]]



class BondAngleComparisonErrorModel(ComparisonErrorModel):
    
    

    def __init__(self, name):
        super().__init__(name)
        self.DataType = namedtuple('DataEntry', ['atom1', 'a1', 'atom2', 'a2', 'atom3', 'a3', 'ref', 'calc', 'error'])
        
    
    def add_value(self, atom1, index1, atom2, index2, atom3, index3, ref_value, calc_value):
        pair = '{0}-{2}-{1}'.format(*list(sorted([atom1, atom3])), atom2)
        error = calc_value - ref_value
        if (pair in self.data):
            self.data[pair] += [[atom1, index1, atom2, index2, atom3, index3, 
            ref_value, calc_value, error]]
        else:
            self.data[pair] = [[atom1, index1, atom2, index2, atom3, index3,
            ref_value, calc_value, error]]  
    

class StructureComparator:
    
    def __init__(self, mof_info_ref, mof_info_calc):
        self.mof_info_ref = mof_info_ref
        self.mof_info_calc = mof_info_calc
        self.error_models = {}

    def compare(self):
        self.compare_bond_length()
        self.compare_bond_angle()

    def compare_bond_length(self):       
        model = BondLengthComparisonErrorModel("bond-length")
        for a1, a2 in self.mof_info_ref.bond_length_indexes:
            atom1 = self.mof_info_ref.struct.sites[a1].specie.symbol
            atom2 = self.mof_info_ref.struct.sites[a2].specie.symbol
            r1 = self.mof_info_ref.struct.get_distance(a1, a2)
            r2 = self.mof_info_calc.struct.get_distance(a1, a2)
            model.add_value(atom1, a1, atom2, a2, r1, r2)
        self.error_models['bond-length'] = model

    def compare_bond_angle(self):       
        model = BondAngleComparisonErrorModel("bond-angle")
        for a1, a2, a3 in self.mof_info_ref.bond_angle_indexes:
            atom1 = self.mof_info_ref.struct.sites[a1].specie.symbol
            atom2 = self.mof_info_ref.struct.sites[a2].specie.symbol
            atom3 = self.mof_info_ref.struct.sites[a3].specie.symbol
            angle1 = self.mof_info_ref.struct.get_angle(a1, a2, a3)
            angle2 = self.mof_info_calc.struct.get_angle(a1, a2, a3)
            model.add_value(atom1, a1, atom2, a2, atom3, a3, angle1, angle2)
        self.error_models['bond-angle'] = model
    
    def as_dict(self):
        err_models = {k: v.as_dict() for k, v in self.error_models.items()}
        d = {
            'error_models': err_models,
            'reference_geom': self.mof_info_ref.struct.as_dict(),
            'calculated_geom': self.mof_info_calc.struct.as_dict()
        }  
        return d;

for frame in get_poscar_frame('XDATCAR'):
    struct = mg.Structure.from_str(frame, fmt='poscar')
    first_mof_info = MOFStructureInfo(struct)
    break

first_struct = first_mof_info.struct


headers = []
for pair_index in first_mof_info.bond_length_indexes:
    header = '{}({})-{}({})'.format(str(first_struct.sites[pair_index[0]].specie), pair_index[0],
                                    str(first_struct.sites[pair_index[1]].specie), pair_index[1])
    headers.append(header)

fout = open ('M-X.dat', 'w') 
print ('#', ', '.join(headers), file=fout)                                

step = 1
last_struct = struct
for frame in get_poscar_frame('XDATCAR'):
    struct = mg.Structure.from_str(frame, fmt='poscar')
    row = []
    for pair_index in first_mof_info.bond_length_indexes:
        key = (pair_index[0], pair_index[1])
        r = struct.get_distance(pair_index[0], pair_index[1])
        row.append(r)
    print (step, ', '.join(map(str, row)), file=fout)
    step+=1
    last_struct = struct

last_mof_info = MOFStructureInfo(last_struct)


#Comparison
comparator = StructureComparator(first_mof_info, last_mof_info)
comparator.compare()

print('')
print ('='*80)
print('{:^80s}'.format('Bond length analysis'))
print('{:35s}{:>15s}{:>15s}{:>15s}'.format('Entry', 'Reference', 'Calculated', 'Error'))
print ('-'*80)

bl_model = comparator.error_models['bond-length']
if (bl_model):
    for values in bl_model.data.values():
        ref_values = []
        calc_values = []
        for atom1, a1, atom2, a2, ref, calc, err in values:
            header = '{}({})-{}({})'.format(atom1, a1, atom2, a2)
            print('{:35s}{:15.8f}{:15.8f}{:15.8f}'.format(header, ref, calc, err))
            ref_values.append(ref)
            calc_values.append(calc)
        print('-'*80)
        print('{:35s}{:15.8f}{:15.8f}'.format('AVG', statistics.mean(ref_values), statistics.mean(calc_values)))
        print('')



print('')
print ('='*80)
print('{:^80s}'.format('Bond angle analysis'))
print('{:35s}{:>15s}{:>15s}{:>15s}'.format('Entry', 'Reference', 'Calculated', 'Error'))
print ('-'*80)

ba_model = comparator.error_models['bond-angle']
if (ba_model):
    for values in ba_model.data.values():
        for atom1, a1, atom2, a2, atom3, a3, ref, calc, err in values:
            header = '{}({})-{}({})-{}({})'.format(atom1, a1, atom2, a2, atom3, a3)
            print('{:35s}{:15.8f}{:15.8f}{:15.8f}'.format(header, ref, calc, err))

print ('='*80)



print('')
print('{:^80s}'.format('Coordination Numbers'))
print('='*80)
print('{:20s}{:>20s}{:>20s}{:>20s}'.format('Metal', 'Init', 'Final', 'Error'))
print ('-'*80)

cn_first = first_mof_info.coordination_number 
cn_last = last_mof_info.coordination_number
changed = False
for i in range(len(cn_first)):
    mi = first_mof_info.metal_indexes[i]
    header = '{}({})'.format(str(first_struct.sites[mi].specie), mi)
    print('{:20s}{:20d}{:20d}{:20d}'.format(header, cn_first[i], cn_last[i], cn_last[i]-cn_first[i]))
    if ((cn_last[i]-cn_first[i])  != 0):
        changed = True
print('-'*80)
if (changed):
    print('WARNING: CN changed')
print('='*80)



print('')


print('{:^80}'.format('Error Statistics'))
print ('='*80)
stat = bl_model.get_stat_value()
print ('{:14s}{:>6s}{:>15s}{:>15s}{:>15s}{:>15s}'.format('Bond', 'N', 'MAX', 'MSE', 'MAD', 'RMS'))
print ('-'*80)
for key, value in stat['detailed'].items():
     print ('{:14s}{:>6d}{:>15.8f}{:>15.8f}{:>15.8f}{:>15.8f}'.format(key,
     value['npoints'], value['max'], value['mse'], value['mad'], value['rms']))
print ('-'*80)
print ('{:14s}{:>6d}{:>15.8f}{:>15.8f}{:>15.8f}{:>15.8f}'.format('All-BL',
     stat['npoints'], stat['max'], stat['mse'], stat['mad'], stat['rms']))
print ('='*80)

stat = ba_model.get_stat_value()
print ('{:14s}{:>6s}{:>15s}{:>15s}{:>15s}{:>15s}'.format('Angle', 'N', 'MAX', 'MSE', 'MAD', 'RMS'))
print ('-'*80)
for key, value in stat['detailed'].items():
     print ('{:14s}{:>6d}{:>15.8f}{:>15.8f}{:>15.8f}{:>15.8f}'.format(key,
     value['npoints'], value['max'], value['mse'], value['mad'], value['rms']))
print ('-'*80)
print ('{:14s}{:>6d}{:>15.8f}{:>15.8f}{:>15.8f}{:>15.8f}'.format('All-BA',
     stat['npoints'], stat['max'], stat['mse'], stat['mad'], stat['rms']))
print ('='*80)




with open('analysis.json', 'w') as f:
    print(json.dumps(comparator.as_dict(),indent=2), file=f)





            





