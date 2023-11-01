import os
import parse
from chemlib import Element
import numpy as np
import pandas as pd
from Bio.PDB import *


""" 
Pipeline:
    remove_duplicate_pdbs.py (optional)
    gencsv.py                           <--
    mkinpt.py
    solvate_all_{linux/windows}.py
    batch_run_analysis.py (controller)
        analyze_wires.py (main process)
        
        
For each PDB file in ./*ases/, this code gathers
the data needed to make a packmol inp file, which is 
placed in ./*ases/{model ID}.csv.

Following running this code, the user runs mkinpt.py
to create ./*ases/{model ID}.inp for each PDB.

PACKMOL is then ran by the user to make
    ./*ases/{model ID}_solvated.pdb 
using 
    ./*ases/{model ID}.pdb,
    ./*ases/{model ID}.inp
"""


# Water molecules per A; calculated using 1g/mL
WATER_PER_A3 = 0.03

# Buffer region for artificial solvation sphere.
# The distance between the furthest enzyme atom
# from the center of the enzyme and the boundary
# of the sphere to fill with water. 
BUFFER = 5

# Gather all PDB files from all ./*ases/ dirs
inpt_dirs = [d for d in os.scandir() if d.name[-4:]=='ases']
dir_contents = [f.path for d in inpt_dirs for f in os.scandir(d)]
pdbs = [f for f in dir_contents if '.pdb' in f]


def get_radius(e):
    ''' Call Element to get the radius of e from the
    name of the element as determined by biopython.
    ''' 
    def radius(e):
        e = Element(e)
        return float(e.AtomicRadius)
    
    try:
        return radius(e)    
    except:
        pass
    
    # e may not be recognized by Element.
    # Change all letters except first to lowercase,
    # e.g. CL --> Cl, and remove '\n'.
    e = str(e[0] + e[-1].lower()).replace('\n','')
    try:
        return radius(e)
    # In this case, e may be a marker in the file
    except:
        return None


# For converting an atom's radius to its volume
def rad_to_vol(r):
    return (4/3) * np.pi * r**3
    

for file in pdbs:
    # Remember, file is actually file.path
    newfile = file.replace('.pdb', '.csv')

    # Uncomment to skip files that have been processed
    if os.path.exists(newfile):
        continue
    
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('id', file)
    residues = [r for r in struct.get_residues()]

    # Gather all atoms' elements and coordinates
    # from the PDB file.
    elements = []
    coords_list = []
    
    for residue in residues:
        atoms = list(residue.get_atoms())
        
        # Water gets a special treatment because molecular
        # and elemental naming conventions vary, and it is 
        # often placed in the model by the originators 
        # without Hydrogen.
        
        # If there is only 1 atom AND the molecule is water,
        # then manually add a water molecule's worth of 
        # atoms to the list of elements.
        is_water = residue.get_resname() in ('HOH', 'SOL')
        if is_water and len(atoms) == 1:
            for _ in range(2): 
                elements.append('H')   
            elements.append('O')
            continue
        
        # Add elements of atoms to list of elements in PDB
        # Add coordinates of atoms to list of coords in PDB
        for atom in atoms:
            elements.append(atom.element)
            coords_list.append(atom.coord)

    radii = [get_radius(e) for e in elements]
    
    # Skip files if any elements in them are unparsable
    # Comment out to ignore them.
    if any(r == None for r in radii):
        continue
    
    volumes = (rad_to_vol(r) for r in radii)
    existing_volume = sum(volumes)

    coords_list = np.array(coords_list)
    # Get coords of enzyme, box extremum
    mins, maxes = zip(*[(min(x), max(x)) for x in coords_list.T])

    # Get a decent solvation sphere radius
    middle = np.array(
                      np.array(mins) +
                     (np.array(maxes) - np.array(mins)) / 2
        )
    deviations = np.linalg.norm(middle - coords_list, axis=1)
    radius = int(max(deviations) + BUFFER)


    #### Create the csv file ####
    
    # Done in try/except to avoid killing the
    # process when one file causes an error. 
    try:
        dir, name = os.path.split(file)
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

        df = pd.DataFrame.from_dict(pdb_info, orient='index')
        df.to_csv(newfile)
        print(f"{name} finished processing")
     
    except Exception as e:
        print()
        print("Encountered error while processing", file)
        print(e)
        print("Continuing")
        print()
        continue
