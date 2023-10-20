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

files = []
for dir in os.scandir():

    if not 'ases' in dir.name:
        continue

    for file in os.scandir(dir):
        if '.csv' in file.name:
            files.append(file)

# Make a dict of dicts.
# Each dict contains the values critical to 
# solvation with PACKMOL and will be used to
# fill a copy of model.inp with data.
model_dict = {}
for file in files:
   
    lines = []
    with open(file,'r') as f:
        for i,line in enumerate(f):
            if i > 0:
                lines.append(line)

    for i,line in enumerate(lines):
        for t in ('[', ']',
                  '(', ')',
                  '\n', "\"",):
            line = line.replace(t,'')
        line = line.split(',')
        lines[i] = line
    
    if len(lines) == 0:
        print(f"File {file.path} is empty.")
        print("Deleting this file.")
        os.remove(file)
        continue
        
    # Directory is the first row and 
    # second column in the csv file.
    dir = lines[0][1]
    del lines[0]
    
    # Dict of csv values 
    d = {}
    for line in lines:
        d[line[0]] = int(float(line[1]))
    
    # Model and Path differ in that the
    # latter contains the file extension.
    model = file.name[:-4]
    d['Model'] = os.path.join(dir, model)
    d['Outfile'] = path + '_solvated'
    d['Path'] = file.path

    model_dict[model] = d


with open('model.inp','r') as f:
    standard_model = [l for l in f]
    standard_model = ''.join(standard_model)


# First column in this matrix are the placeholder
# values in standard_model for the actual values in 
# ./*ases/*.csv
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
