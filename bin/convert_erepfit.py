#!/usr/bin/env python3

import json
import yaml

import sys

input = sys.argv[1]

with open(input, 'r') as f:  
    doc = json.load(f)

del doc["version"]
del doc['options']['use_total_minus_repulsive']


for i, system in enumerate(doc['systems']):
    shifts = system['geometry']['kpoint_shifts']
    kpoints = system['geometry']['kpoints']
    del system['geometry']['kpoint_shifts']
    del system['geometry']['kpoints']
    new_kpoints = {
        'supercell_folding': [
            [kpoints[0], 0, 0],
            [0, kpoints[1], 0],
            [0, 0, kpoints[2]],
            shifts
        ]
    }
    system['geometry']['kpoints'] = new_kpoints

with open('erepfit_input_new.json', 'w') as fout:
    print(json.dumps(doc, indent=2, sort_keys=True), file=fout)

with open('erepfit_input_new.yaml', 'w') as fout:
    yaml.dump(doc, fout)