from parse import *
import os
import numpy as np
from itertools import product
import random
import argparse
import time
import pandas as pd
from Bio.PDB import *


''' This program takes pdb files {ID}_solvated.pdb that have been solvated
by packmol, extracts active site locations from corresponding .index files
which contain active site line numbers in the ORIGINAL pdb file {ID}.pdb,
writes the active site lines into their own pdb, then reads this pdb to
obtain active site residues in the solvated file. Both the solvated file and
the temporary active site pdb file are parsed with Bio.PDB.PDBParser. The
program then establishes water in the solvated pdb as SOLVENT or NEAR ENZYME.
It then forms a wire of waters within distance CUTOFF of the active site
atoms and propagates this wire, finding more and more waters within CUTOFF
of each other. It then asks how many solvent water molecules are within
CUTOFF of any component of the wire. This and more data are written to a
file in the target directory, {ID}_wireInfo_cutoff_{CUTOFF}.csv.
'''

p = argparse.ArgumentParser()
p.add_argument('--dir')
p.add_argument('--cutoff',type=float)
args = p.parse_args()
dir = args.dir
CUTOFF = float(args.cutoff)
parser = PDBParser(QUIET=True)


def dist(a,b):
    ''' Euclidian distance between cartesian coordinates.
    '''
    r = 0
    for i in range(3):
        r += (a[i] - b[i]) ** 2
    return r**0.5


def relatively_close(X, Y, cutoff):
    ''' Exists to speed up nested looping for assessing
    whether distance between X and Y are less than
    cutoff.
    '''
    x, y, z = X
    i, j, k = Y
    if abs(x - i) > cutoff:
        return False
    if abs(y - j) > cutoff:
        return False
    if abs(z - k) > cutoff:
        return False
    else:
        return True


def filter_by_distance(arr_1, arr_2, cutoff):
    ''' Return arr_1 values that are within cutoff
    of ANY arr_2 values.
    '''
    arr_1_within_distance = []
    for X in arr_1:
        for Y in arr_2:
            if not relatively_close(X, Y, cutoff):
                continue
            elif dist(X, Y) < cutoff:
                arr_1_within_distance.append(X)
                break
    return arr_1_within_distance


class PDB_Object():
    def __init__(self,file):
        self.file = file
        self.id = file.name[:-4]
        self.solvated_file = os.path.join(dir, self.id + '_solvated.pdb')
        self.index_file = os.path.join(dir, self.id + '.index')
        self.active_site_pdb = self.id + '_active_site.pdb'
        assert os.path.exists(self.solvated_file)

    def get_as_line_numbers(self):
        ''' Get lines numbers of active site atoms from
        index file. Returns after first line because index
        file only has 1 line.
        '''
        with open(self.index_file,'r') as f:
            for line in f:
                self.as_line_numbers = str2list(line)
                return

    def has_index_file(self):
        ''' Does the file containing AS line numbers
        exists and is it not empty?
        '''
        if os.path.exists(self.index_file):
            self.get_as_line_numbers()
            if len(self.as_line_numbers) != 0:
                return True
            else:
                return False
        else:
            return False

    def get_active_site_lines(self):
        ''' Extract lines from pdb file containing the active
        site atoms whose postions are given in index file.
        '''
        with open(self.file, 'r') as f:
            enum = enumerate([line for line in f])
            self.as_lines = [l for i, l in enum if i in self.as_line_numbers]

    def write_active_site_file(self):
        ''' Write file containing just active site lines.
        '''
        with open(pdb.active_site_pdb, 'w') as f:
            for line in self.as_lines:
                f.write(line)

    def bioparse_file(self,file):
        ''' Return residues (bio module) in file.
        '''
        struct = parser.get_structure(self.id, file)
        return [r for r in struct.get_residues()]

    def remove_active_site_file(self):
        ''' Delete active site file for space.
        '''
        os.system('rm {}'.format(self.active_site_pdb))

    def get_active_site_residues(self):
        ''' Get active site residues from active site file.
        '''
        self.act_st_resids = self.bioparse_file(self.active_site_pdb)

    def get_solvated_file_residues(self):
        ''' Get residues from file solvated by packmol.
        '''
        self.residues = self.bioparse_file(self.solvated_file)

    def get_outfile_path(self):
        ''' Return path of output file to make.
        '''
        outfile_template = '{}_wireInfo_cutoff_{}.csv'
        outfile = outfile_template.format(self.id, int(CUTOFF))
        return os.path.join(dir, outfile)

    def has_outfile(self):
        ''' Check if outfile exists.
        '''
        outfile_path = self.get_outfile_path()
        return os.path.exists(outfile_path)


def pdb_not_solvated(file):
    ''' Want pdb file, but not solvated pdb file.
    The former will be used to ID active site residues.
    '''
    return 'pdb' in file.name and not 'solvated' in file.name


''' Get viable pdb files for testing. These have versions
solvated by packmol (checked by assertion in class), index files
containing the active site, and NO outfiles previously made by this
program. '''
pdbs_in_dir = [f for f in os.scandir(dir) if pdb_not_solvated(f)]
pdb_file_candidates = [PDB_Object(f) for f in pdbs_in_dir]
pdb_files = [f for f in pdb_file_candidates if f.has_index_file()]
untested_pdb_files = [f for f in pdb_files if not f.has_outfile()]

# Main operation on each pdb file
for pdb in untested_pdb_files:
    pdb.get_active_site_lines() # get AS from orig PDB
    pdb.write_active_site_file() # make file of just AS
    pdb.get_active_site_residues() # read file of just AS
    pdb.remove_active_site_file() # remove file of just AS
    pdb.get_solvated_file_residues() # read solvated file

    # sort residues into active site, water, enzyme (not active site)
    active_site_coords, water, enz_coords = [], [], []
    for r in pdb.residues:
        for atom in r.get_atoms():
            c = list(atom.coord)
            if r.resname in ('HOH','SOL') and atom.element=='O':
                water.append(c)
            elif r in pdb.act_st_resids:
                active_site_coords.append(c)
            else:
                enz_coords.append(c)

    # get water near enzyme
    hoh_near_enz = filter_by_distance(water, enz_coords, CUTOFF)
    # water not near enzyme
    solvent = [h for h in water if not h in hoh_near_enz]
    # seed wire -- water near active site
    wire = filter_by_distance(hoh_near_enz, active_site_coords, CUTOFF)
    # make water near enz and nascent wire mutually exclusive
    hoh_near_enz = [h for h in hoh_near_enz if not h in wire]

    # Expand wire by CUTOFF until it cannot be further expanded
    while True:
        wire_mask = wire
        for x in wire:
            for i, y in enumerate(hoh_near_enz):
                if not relatively_close(x, y, CUTOFF):
                    continue
                elif dist(x, y) < CUTOFF:
                    del hoh_near_enz[i]
                    wire.append(y)
        # if we haven't added anything new to wire
        if wire_mask is wire:
            break

    solv_contacted = filter_by_distance(solvent, wire, CUTOFF)

    keys = ['Enzyme Atoms',
            'AS Atoms',
            'Wire Length',
            'Water Near enz_coords',
            'N Solvent',
            'N Solvent Contacted']

    values = [len(enz_coords),
              len(active_site_coords),
              len(wire),
              len(hoh_near_enz) + len(wire),
              len(solvent),
              len(solv_contacted)]

    df = pd.DataFrame(keys, values)
    df.to_csv(pdb.get_outfile_path())
