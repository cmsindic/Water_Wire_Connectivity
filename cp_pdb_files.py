import os
from shutil import copyfile

dirs = [d for d in os.scandir('../') if 'ases' in d.name[-4:]]
for dir in dirs:
    new_dir = dir.name.replace('..','')
    os.mkdir(new_dir)
files = [obj for d in dirs for obj in os.scandir(d) if 'pdb' in obj.name]
for pdb in files:
    outfile = pdb.path.replace('..','.')
    copyfile(pdb.path, outfile)
