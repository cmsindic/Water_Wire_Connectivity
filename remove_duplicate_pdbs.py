import os
import re
from itertools import *


""" 
Pipeline:
    remove_duplicate_pdbs.py (optional) <--
    gencsv.py                           
    mkinpt.py
    solvate_all_{linux/windows}.py
    batch_run_analysis.py (controller)
        analyze_wires.py (main process)
        
Remove duplicates from the pool of PDB files
based on the compound name in each PDB.

It is assumed that the PDB files are in 
directories located in the working directory
and that they are grouped by enzyme class, 
meaning that the names of the directories
contain "ases".
"""


# Gather all PDB files from all ./*ases/ dirs
inpt_dirs = [d for d in os.scandir() if d.name[-4:]=='ases']
dir_contents = [f.path for d in inpt_dirs for f in os.scandir(d)]
pdbs = [f for f in dir_contents if '.pdb' in f]

name_dict = {pdb: " " for pdb in pdbs}
for file in pdbs:
    # Parse PDB as a list of lines
    # ... but only lines containing "COMPND"
    # and "ASE"
    with open(file, 'r') as f:
        lines = [l for l in f if 'COMPND' in l]
        lines = [l.replace('\n', '') for l in lines]
        lines = [l.replace(';', '') for l in lines]
        lines = [l for l in lines if 'ASE ' in l]

    # Skip files with no compound name
    if len(lines) == 0:
        continue 
        
    # The string between linestart and ":" 
    # contains the name of the enzyme. This
    # is the same for every file.
    nameline = lines[0]
    to_rm = re.search('^(.*)\:', name_line)
    if to_rm is None:
        continue
    to_rm = to_rm.group(0)

    # May be removing space between linestart
    # and ":" IN THE DESCRIPTOR instead of that
    # between linestart and first occurence of ":"
    # outside of the descriptor.
    # This fixes this issue
    further_remove = re.search('^(.*)\:', to_rm)
    if further_remove is None:
        pass
    else:
        to_rm = further_remove.group(0)
    
    # Remove everything around the possible 
    # enzyme compound name. 
    name_line = name_line.replace(to_rm, '')
    name_line = name_line.split(' ')
    name_line = [n for n in name_line if n != '']
    name_dict[file] = name_line


# Remove all nameless PDBs
name_dict = {k: v for k, v in name_dict.items() if v != ' '}

# Convert the dict to a list so it can sort
nd = [(k, v) for k, v in name_dict.items()]
nd.sort(key=lambda x: x[1])

# Delete all but the first of files that 
# have identical enzyme names.
for k,v in groupby(nd, key=lambda x: x[1]):
    copies = [x[0] for x in v]
    if len(copies) > 1:
        for f in copies[1:]:
            os.remove(f)

