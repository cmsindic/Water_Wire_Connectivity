import os
import re
from itertools import *

classDirs = [d for d in os.scandir('..') if d.name[-4:]=='ases']
#classDirs = [d for d in classDirs if 'hydrolases' in d.name]

pdbs = []
for dir in classDirs:
    try:
        os.mkdir(dir.name)
    except:
        pass
    for file in os.scandir(dir):
        if '.pdb' in file.name:
            pdbs.append(file.path)

name_dict = {pdb:" " for pdb in pdbs}

for file in pdbs:

    with open(file, 'r') as f:
        lines = [l for l in f if 'COMPND' in l]

    lines = [l.replace('\n', '') for l in lines]
    lines = [l.replace(';', '') for l in lines]
    lines = [l for l in lines if 'ASE ' in l]

    try:
        name_line = lines[0]
    except IndexError:
        continue

    # string between start and : -- same for every filecombinations
    to_rm = re.search('^(.*)\:', name_line)
    if to_rm is None:
        continue
    to_rm = to_rm.group(0)

    # may be removing space between start and : IN DESCRIPTOR
    # instead of that betwen start and first occurence of :
    # this fixes this issue
    further_remove = re.search('^(.*)\:', to_rm)
    if further_remove is None:
        pass
    else:
        to_rm = further_remove.group(0)

    name_line = name_line.replace(to_rm, '')
    name_line = name_line.split(' ')
    name_line = [n for n in name_line if n!='']
    name_dict[file] = name_line

name_dict = {k:v for k,v in name_dict.items() if v != ' '}

sames = []
for k,v in name_dict.items():
    for kk, vv in name_dict.items():
        if k == kk:
            continue
        if v == vv:
            sames.append([k,kk,v])

nd = [(k,v) for k,v in name_dict.items()]
nd.sort(key=lambda x: x[1])

for k,v in groupby(nd, key=lambda x: x[1]):
    copies = [x[0] for x in v]
    if len(copies) > 1:
        for f in copies[1:]:
            os.remove(f)

'''
files = [x[0] for x in v]
if len(files) > 1:
    for file in files[1:]:
        os.remove(file)'''
