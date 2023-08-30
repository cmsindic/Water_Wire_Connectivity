import os
import pandas as pd
import numpy as np

'''collect data from csv files made by gencsv.py and
use their data to make an inp for packmol for each pdb
model using model.inp as a framework'''

tr = ['[',']','(',')','\n',"\""]

files = []
for dir in os.scandir():

    if not 'ases' in dir.name:
        continue

    for file in os.scandir(dir):
        if '.csv' in file.name:
            files.append(file)

modelDict = {}
for file in files:

    model = file.name[:-4]

    lines = []
    with open(file,'r') as f:
        for i,line in enumerate(f):
            if i > 0:
                lines.append(line)

    for i,line in enumerate(lines):
        for t in tr:
            line = line.replace(t,'')
        line = line.split(',')
        lines[i] = line

    dir = lines[0][1]
    del lines[0]

    dict = {}
    for line in lines:
        dict[line[0]] = int(float(line[1]))

    path = os.path.join(dir,model)
    dict['Model'] = '../' + path
    dict['Outfile'] = path + '_solvated'
    dict['Path'] = file.path

    modelDict[model] = dict

stdinModel = []
with open('model.inp','r') as f:
    for line in f:
        stdinModel.append(line)
stdinModel = ''.join(stdinModel)

replaceList = [
    ['__MODEL__','Model'],
    ['__NWAT__','Number of Waters'],
    ['__RADIUS__','Radius'],
    ['__OUTFILE__','Outfile']
]

for m,v in modelDict.items():
    inp = stdinModel
    for x in replaceList:
        inp = inp.replace(x[0],str(v[x[1]]))
    inpName = v['Path'][:-4] + '.inp'
    print(inpName)
    with open(inpName,'w') as f:
        f.write(inp)
