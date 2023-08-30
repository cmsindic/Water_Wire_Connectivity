import os
import pandas as pd

fname_key = 'wireInfo_cutoff5_new'
dirs = [d for d in os.scandir() if 'ases' in d.name]
contents = [f for dir in dirs for f in os.scandir(dir)]
data_files = [f for f in contents if fname_key in f.name]

mega = {'Enzyme Atoms':[],
        'AS Atoms': [],
        'Wire Length': [],
        'Water Near enz_coords': [],
        'N Solvent': [],
        'N Solvent Contacted': [],
        'Class': [],
        'Filename': []}

for file in data_files:

    with open(file, 'r') as f:
        lines = [l for l in f]

    if len(lines) == 0: continue

    del lines[0]

    lines = [l.replace('\n','') for l in lines]
    lines = [l.split(',') for l in lines]

    for value, key in lines:
        mega[key].append(int(value))

    dir, name = os.path.split(file.path)
    name = name[:4]
    mega['Class'].append(dir)
    mega['Filename'].append(name)

df = pd.DataFrame(mega)
df.to_csv("main_out_cutoff5.csv")
