import os
import parse
from chemlib import Element
import numpy as np
import pandas as pd
from Bio.PDB import *

"""for each pdb file in ../*ases/, gather data needed
to make a packmol inp file, and place in ./*ases/*csv
for creation of XXX.inp by mkinpt.py"""

WATER_PER_A3 = 0.03
BUFFER = 5

classDirs = [d for d in os.scandir('..') if d.name[-4:]=='ases']
#classDirs = [d for d in classDirs if not 'ligases' in d.name]
classDirs = [d for d in classDirs if 'oxidoreductases' in d.name]

pdbs = []
for dir in classDirs:
    try:
        os.mkdir(dir.name)
    except:
        pass
    for file in os.scandir(dir):
        if '.pdb' in file.name:
            pdbs.append(file.path)

def rad_to_vol(r):
    return (4/3) * np.pi * r**3

def get_radius(e):
    try:
        return float(Element(e).AtomicRadius)
    except:
        # CL --> Cl & rm \n
        e = str(e[0] + e[-1].lower()).replace('\n','')
        try:
            return float(Element(e).AtomicRadius)
        except:
            return None

for file in pdbs:

    newpath = file.replace('../','') # rm ../ in beginning of path
    dir, name = os.path.split(newpath)
    newfile = newpath[:-4] + '.csv'

    if os.path.exists(newfile):
        continue

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('id', file)
    residues = [r for r in struct.get_residues()]

    elements = []
    for residue in residues:
        atoms = list(residue.get_atoms())
        if residue.get_resname() in ('HOH', 'SOL'):
            if len(atoms) == 1:
                for _ in range(2): elements.append('H')
                elements.append('O')
                continue
        for atom in atoms:
            elements.append(atom.element)

    radii = [get_radius(e) for e in elements]
    if any(r == None for r in radii):
        continue
    existing_volume = sum([rad_to_vol(r) for r in radii])
    coords_list = [a.coord for r in residues for a in r.get_atoms()]

    try:


        coords_list = np.array(coords_list)
        # get coords of enzyme, box extremum
        mins, maxes = zip(*[(min(x), max(x)) for x in coords_list.T])

        # get a decent solvation sphere radius
        middle = np.array(
                          np.array(mins) +
                         (np.array(maxes) - np.array(mins)) / 2
            )
        deviations = np.linalg.norm(middle - coords_list, axis=1)
        radius = int(max(deviations) + BUFFER)

        enz_volume = round(existing_volume, 2)
        total_volume = (4/3) * np.pi * radius**3
        water_volume = round(total_volume - enz_volume, 2)
        number_of_water = int(water_volume * WATER_PER_A3)

        pdb_info = {
            'Directory':dir,
            'Enzyme Volume':enz_volume,
            'Total Volume':total_volume,
            'Water Volume':water_volume,
            'Number of Waters':number_of_water,
            'Radius':radius
        }

        #print(pdb_info)

        df = pd.DataFrame.from_dict(pdb_info,orient='index')
        df.to_csv(newfile)
        print(name)
    except:
        print('No COORDS FOR', file)
