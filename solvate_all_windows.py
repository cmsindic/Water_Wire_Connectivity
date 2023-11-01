import os

""" 
Pipeline:
    remove_duplicate_pdbs.py (optional) 
    gencsv.py                           
    mkinpt.py                           
    solvate_all_{linux/windows}.py     <--
    batch_run_analysis.py (controller)
        analyze_wires.py (main process)

Runs packmol on all inp/pdb file pairs in ./*ases
that have not already been solvated. 

Works for windows. 
"""

dirs = [d for d in os.scandir() if 'ases' in d.name]
models = set((d.name, f.name[:4]) for d in dirs for f in os.scandir(d))

for d, i in models:

    sfile = i +'_solvated.pdb'
    sf_path = os.path.join(d, sfile)
    if os.path.exists(sf_path):
        continue

    inpt_file = i + '.inp'
    inpt_path = os.path.join(d, inpt_file)
    if not os.path.exists(inpt_path):
        continue
    
    pdb = i + '.pdb'
    pdb_path = os.path.join(d, pdb)
    cmd = "packmol < {}".format(inpt_path)
    os.system(cmd)
