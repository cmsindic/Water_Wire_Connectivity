import os
from shutil import copy


dirs = [d for d in os.scandir('backup') if 'ases' in d.name]
files = [f for d in dirs for f in os.scandir(d) if 'no_info' in f.name]

for f in files:
    op = f.path
    np = op.replace('backup','.')
    open(np, 'w').close()
    copy(op, np)
