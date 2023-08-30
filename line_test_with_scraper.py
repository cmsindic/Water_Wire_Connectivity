import parse
import os
import argparse
import pandas as pd
from Bio.PDB import *
from lxml import etree
import requests
from io import StringIO
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from time import time
import numpy as np



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


class PDB_Object():
    def __init__(self,file):
        self.file = file
        self.name = file.name
        name, ext = file.name.split('.')
        self.id = name[:4]
        self.no_as_file = dir + '/' + self.id + '.no_info'
        self.no_solv_file = dir + '/' + self.id + '.no_solv'

    def get_residues(self):
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(self.id, self.file)
        self.residues = [r for r in struct.get_residues()]

    def get_rcsb_tree(self):
        ''' Return the xml tree of the rscb page.
        '''
        parser = etree.HTMLParser()
        base_rcsb_url = 'https://www.rcsb.org/structure/'
        url = base_rcsb_url + self.id
        print(url)
        page = requests.get(url)
        html = page.content.decode("utf-8")
        return etree.parse(StringIO(html), parser=parser)

    def check_for_uniprot_id(self):
        ''' Resolve the identity of the uniprot object linked
        to the pdb model by extracting the name from the uniprot
        url that is linked in the rscb page for the model.
        '''
        tree = self.get_rcsb_tree()
        refs = tree.xpath("//a")
        links = [link.get('href', '') for link in refs]
        uniprot_links = [l for l in links if 'uniprot.org' in l]
        if len(uniprot_links) == 0:
            return False
        # Ex: 'https://www.uniprot.org/uniprot/P00942' --> 'P00942'
        self.uniprot_id = uniprot_links[0].split('/')[-1]
        return True

    def fetch_uniprot_text(self):
        ''' Get the raw text of the uniprot text file corresponding
        with the pdb model.
        '''
        uniprot_base_url = 'https://rest.uniprot.org/uniprotkb/'
        uniprot_id = self.uniprot_id
        url = uniprot_base_url + uniprot_id + '.txt'
        page = requests.get(url)
        return page.text

    def resolve_res_ids(self):
        ''' Extract the residue ID numbers of the active site
        residues from the uniprot text file for the pdb model.
        '''
        # Text file from uniprot detailing model features
        raw_text = self.fetch_uniprot_text()
        # Split text into list of lines
        text = raw_text.split('\n')
        # Split each line into list of words
        parsed = (parse.text_to_list(line) for line in text)
        # Isolate lines containing active site info
        act_site_lines = (line for line in parsed if 'ACT_SITE' in line)
        # Residue IDs are the last words in the lines
        return tuple(int(line[-1]) for line in act_site_lines)

    def get_active_site_residues(self):
        ''' Get residues belonging to the active site.
        '''
        act_site_ids = self.resolve_res_ids()
        self.as_res = [r for r in self.residues if r.id[1] in act_site_ids]
        self.as_res = [r for r in self.as_res if not r.resname in ('HOH','SOL')]

    def has_active_site_information(self):
        return len(self.as_res) > 0

    def get_outfile_path(self):
        ''' Return path of output file to make.
        '''
        outfile_template = '{}_wireInfo_cutoff{}_new.csv'
        outfile = outfile_template.format(self.id, int(CUTOFF))
        return os.path.join(dir, outfile)

    def make_no_as_file(self):
        ''' Make a file to indicate that the model has no
        active site information'''
        os.system('touch {}'.format(self.no_as_file))

    def make_no_solv_file(self):
        ''' Make a file to indicate that the model has no
        active site information'''
        os.system('touch {}'.format(self.no_solv_file))

    def has_as_file(self):
        ''' Check if file exists to indicate that no AS info can
        be found.
        '''
        return os.path.exists(self.no_as_file)


def members_in_proximity(p1, p2):
    ''' Get elements of p1 within CUTOFF of
    any members of p2.
    '''
    MAX_MAT_SIZE = 5 * 10**8
    mat_size = len(p1) * len(p2)

    CHUNK_SIZE = 2 * 10**3
    n_full_chunks = int(len(p1) / CHUNK_SIZE) + 1
    chk_inds = [c * CHUNK_SIZE for c in range(n_full_chunks)] + [None]
    iter_chnks = range(len(chk_inds) - 1)
    split_arr = [p1[chk_inds[i]:chk_inds[i + 1]] for i in iter_chnks]

    if mat_size < MAX_MAT_SIZE:
        d = enumerate(distance_matrix(p1, p2))
        return [p1[i] for i,r in d if any(r < CUTOFF)]

    else:
        print(mat_size)
        inds = []
        for ii,p in enumerate(split_arr):
            d = enumerate(distance_matrix(p, p2))
            inds += [i+ii for i,r in d if any(r < CUTOFF)]
        return [p1[iii] for iii in inds]


def pdb_solvated(file):
    ''' Want pdb file, but not solvated pdb file.
    The former will be used to ID active site residues.
    '''
    return 'pdb' in file.name and 'solvated' in file.name


# Get viable pdb files for testing. These are solvated by packmol.
pdbs_in_dir = [f for f in os.scandir(dir) if pdb_solvated(f)]
pdb_files = [PDB_Object(f) for f in pdbs_in_dir]
pdb_files = [f for f in pdb_files if not f.has_as_file()]
pdb_files = [f for f in pdb_files if not os.path.exists(f.get_outfile_path())]

# Main operation on each pdb file
for pdb in pdb_files:
    t1 = time()
    print(pdb.id)

    if not pdb.check_for_uniprot_id():
        pdb.make_no_as_file()
        continue

    pdb.get_residues()
    pdb.get_active_site_residues()

    if len(pdb.as_res) == 0:
        pdb.make_no_as_file()
        continue

    print("Active site residues:", pdb.as_res)

    # sort residues into active site, water, enzyme (not active site)
    active_site_coords, water, enz_coords = [], [], []
    for r in pdb.residues:
        for atom in r.get_atoms():
            c = list(atom.coord)
            if r.resname in ('HOH','SOL') and atom.element=='O':
                water.append(c)
            elif r in pdb.as_res:
                active_site_coords.append(c)
            elif r.resname not in ('HOH','SOL'):
                enz_coords.append(c)

    # get water near enzyme
    hoh_near_enz = members_in_proximity(water, enz_coords)
    # water not near enzyme
    solvent = [h for h in water if not h in hoh_near_enz]
    # seed wire -- water near active site
    # get water near enzyme
    wire = members_in_proximity(hoh_near_enz, active_site_coords)
    # make water near enz and nascent wire mutually exclusive
    n_h_near_enz = len(hoh_near_enz)

    if len(wire) !=0:
        # Expand wire by CUTOFF until it cannot be further expanded
        wire_mask = wire
        while True:
            hoh_near_enz = [h for h in hoh_near_enz if not h in wire]
            if len(hoh_near_enz) == 0: break
            wire_mask = members_in_proximity(hoh_near_enz, wire_mask)
            wire += wire_mask
            if len(wire_mask) == 0: break

    if len(wire) != 0 and len(solvent) != 0:
        solv_contacted = members_in_proximity(solvent, wire)
    else:
        solv_contacted = []

    keys = ['Enzyme Atoms',
            'AS Atoms',
            'Wire Length',
            'Water Near enz_coords',
            'N Solvent',
            'N Solvent Contacted']

    values = [len(enz_coords),
              len(active_site_coords),
              len(wire),
              n_h_near_enz,
              len(solvent),
              len(solv_contacted)]

    print({keys[i]:values[i] for i in range(len(keys))})
    df = pd.DataFrame(keys, values)
    df.to_csv(pdb.get_outfile_path())
    print(time()-t1)

    del pdb
