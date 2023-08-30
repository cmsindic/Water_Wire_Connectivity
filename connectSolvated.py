from parse import *
import os
import numpy as np
from itertools import product
import random
import argparse
import time
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--dir')
#p.add_argument('--cutoff',type=float)
args = p.parse_args()
dir = args.dir
#cutoff = float(args.cutoff)

def dist(a,b):
    r = 0
    for i in range(3):
        r += ( a[i] - b[i] ) ** 2
    r = np.sqrt(r)
    return r

for file in os.scandir(dir):

    # PDB ID
    id = file.name[:-4]

    pdb = os.path.join(dir,id + '_solvated' + '.pdb')
    indexFile = os.path.join(dir,id + '.index')

    solvated = os.path.exists(pdb)
    hasIndex = os.path.exists(indexFile)

    if not (solvated and hasIndex):
        continue

    outfile = os.path.join(dir,id+'_wireInfo_cutoff5.csv')
    if os.path.exists(outfile):
        continue

    # get line numbers of AS atoms
    lineNumbers = []
    with open(indexFile,'r') as f:
        for line in f:
            lineNumbers = str2list(line)

    # if no active site line numbers
    # OR no index file (bc lineNumbers is init as [])
    if len(lineNumbers) == 0:
        continue

    # need only be immutable now
    lineNumbers = tuple(lineNumbers)

    # init 1x3 arrays
    asc,wat,enz = [],[],[]

    anisou = False
    # get enz of AS atom from line numbers
    with open(pdb,'r') as f:
        for i,line in enumerate(f):

            if anisou:
                continue

            # don't process unnecessary lines
            if not ('ATOM' in line or 'HETATM' in line):
                continue

            line = pdbline(line)

            if 'H' in line:
                #print(line)
                continue

            #if 'HW'

            # only want 'neat' lines
            if len(line) != 12:
                if 'HOH' in line:
                    if len(line) == 9 or len(line) == 8:
                        c = tuple(float(x) for x in line[-3:])
                    else:
                        print(line)
                        continue
                elif len(line) == 13:
                    line[2] = line[2] + line[3]
                    del line[3]
                    c = tuple(float(x) for x in line[6:9])

            else:
                try:
                    c = tuple(float(x) for x in line[6:9])
                except:
                    continue

            # if enzyme atom line
            if 'ATOM' in line[0]:

                # add coords to list of enzyme lines
                enz.append(c)

                # add AS coords to special array
                if i in lineNumbers:

                    if 'ANISOU' in line:
                        anisou = True
                    else:
                        asc.append(c)

            # if water atom line
            elif ('HOH' in line[3]):
                wat.append(c)

    if anisou:
        continue

    if len(asc) == 0:
        continue

    solvent = []
    waterNearEnzyme = []
    for w in wat:
        near = 0
        for e in enz:
            if near > 0:
                continue
            if -5 < e[0]-w[0] < 5:
                if -5 < e[1]-w[1] < 5:
                    if -5 < e[2]-w[2] < 5:
                        if dist(w,e) < 5:
                            near += 1
        if near < 1:
            solvent.append(w)
        else:
            waterNearEnzyme.append(w)


    if len(waterNearEnzyme) == 0:
        continue

    wire = []
    indices = []
    waterNearEnzymeMask = waterNearEnzyme
    for i,w in enumerate(waterNearEnzymeMask):
        for e in asc:
            if dist(w,e) < 5:
                wire.append(w)
                indices.append(i)

    waterNearEnzyme = [waterNearEnzyme[i] for i in range(len(waterNearEnzyme)) if not i in indices]
    #wire1 = [x for x in wire]

    # does wire reach bulk?
    while True:
        wire2 = wire
        for x in wire:
            for i,y in enumerate(waterNearEnzyme):
                if dist(x,y) < 5:

                    # unique elements only
                    if not y in wire:

                        # wire and wat are mutually exclusive
                        del waterNearEnzyme[i]
                        wire.append(y)


        # if we haven't added anything new to wire
        if wire2 is wire:
            break

    contacted = set()
    for w in wire:
        for i,s in enumerate(solvent):
            if dist(w,s) < 5:
                contacted.update([i])
    contactedSolvent = [solvent[i] for i in contacted]

    contacted = len(contacted)

    nenz = len(enz)
    nasc = len(asc)
    nwire = len(wire)
    nwat = len(waterNearEnzyme) + len(wire)
    nsolv = len(solvent)
    print(id)
    df = pd.DataFrame(['Enzyme Atoms','AS Atoms','Wire Length','Water Near Enz','N Solvent','N Solvent Contacted'],[nenz,nasc,nwire,nwat,nsolv,contacted])
    df.to_csv(outfile)

    '''
    for i,x in enumerate(wire):
        for y in wire1:
            if x == y:
                del wire[i]

    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    enz = np.array(enz)
    asc = np.array(asc)
    contactedSolvent = np.array(contactedSolvent)
    solvent = np.array(solvent)
    wire = np.array(wire)
    waterNearEnzyme = np.array(waterNearEnzyme)
    ax.scatter3D(asc[:,0], asc[:,1], asc[:,2],
                 color='red', alpha=0.3, s=200)
    ax.scatter3D(enz[:,0], enz[:,1], enz[:,2],
                 color='grey', s=200,alpha=0.01)
    ax.scatter3D(wire[:,0],wire[:,1],wire[:,2],color='green',)#alpha=0.2)
    ax.scatter3D(contactedSolvent[:,0],contactedSolvent[:,1],contactedSolvent[:,2],color='blue',alpha=0.02)
    #fig.savefig('figs/'+args.dir+'_'+id+'.png')
    plt.show()
    '''
