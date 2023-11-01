import os
import pandas as pd
import numpy as np

""" 
Pipeline:
    remove_duplicate_pdbs.py (optional) 
    gencsv.py                           
    mkinpt.py                           <--
    solvate_all_{linux/windows}.py
    batch_run_analysis.py (controller)
        analyze_wires.py (main process)
        
This code takes csv files made by gencsv.py
and converts them into inp files for 
PACKMOL to use to solvate the PDB files 
corresponding with each csv file. 

It requires model.inp to be in the working
directory to be used as a template for its 
output.
"""

def process(line):
    ''' Convert a csv row, as a string,
    to a list of values.
    '''
    fluff = ('[', ']', '(', ')', '\n', '\"')
    for char in fluff:
        line = line.replace(char, '')
    return line.split(',')

# Get all csv files in directories containing PDBs
# Assumes PDBs are stored in directories by class
# of the enzymes, and thus has names *ases.
dirs = [d for d in os.scandir() if "ases" in d.name]
all_files = [f for d in dirs for f in os.scandir(d)]
files = [f for f in all_files if ".csv" in f.name]

# Make a dict of dicts.
# Each dict contains the values critical to 
# solvation with PACKMOL and will be used to
# fill a copy of model.inp with data.
model_dict = {}
for file in files:
    
    # Get lines from csv
    with open(file,'r') as f:
        # Ignore header row in csv 
        header, *lines = [process(l) for l in f]
    
    if len(lines) == 0:
        print(f"File {file.path} is empty.")
        print("Deleting this file.")
        os.remove(file)
        continue
        
    # Directory is the first row and 
    # second column in the csv file.
    direct = lines[0][1]
    del lines[0]
    
    # Dict of csv values 
    d = {l[0]: int(float(l[1])) for l in lines}
    
    # Model and path differ in that the
    # latter contains the file extension.
    model = file.name[:-4]
    path = os.path.join(direct, model)
    d['Model'] = os.path.join(direct, model)
    d['Outfile'] = path + '_solvated'
    d['Path'] = file.path
    
    # Add the dict to the dict of models' vals
    # indexed by the model name of each.
    model_dict[model] = d


with open('model.inp', 'r') as f:
    standard_model = [l for l in f]
    standard_model = ''.join(standard_model)


# First column in this matrix contains the 
# placeholder values in standard_model (model.inp)
# for the actual values in ./*ases/*.csv.
# The second column contains the keys to the dict
# containing those values as they were parsed from 
# each csv file.
replace_list = [
    ['__MODEL__','Model'],
    ['__NWAT__','Number of Waters'],
    ['__RADIUS__','Radius'],
    ['__OUTFILE__','Outfile']
]

# Make a new version of the standard_model
# and edit it with the new data from each csv
# file before saving into ./*ases/.
for v in model_dict.values():
    inp = standard_model
    for x in replace_list:
        inp = inp.replace(x[0], str(v[x[1]]))
        
    input_name = v['Path'][:-4] + '.inp'
    with open(input_name,'w') as f:
        f.write(inp)
    
    print(f"{input_name} created")
