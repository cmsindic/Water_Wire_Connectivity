import os

""" 
Pipeline:
    remove_duplicate_pdbs.py (optional) 
    gencsv.py                           
    mkinpt.py                           
    solvate_all_{linux/windows}.py     <--
    batch_run_analysis.py (controller)
        analyze_wires.py (main process)

Runs packmol on all files in ./*ases that have
not already been solvated. Works for linux. """

dirs = [d for d in os.scandir() if 'ases' in d.name]
models = set((d.name, f.name[:4]) for d in dirs for f in os.scandir(d))

for dir, id in models:

    sfile = id +'_solvated.pdb'
    sf_path = os.path.join(dir, sfile)
    if os.path.exists(sf_path):
        continue

    inpt_file = id + '.inp'
    inpt_path = os.path.join(dir, inpt_file)
    if not os.path.exists(inpt_path):
        continue

    pdb = id + '.pdb'
    pdb_path = os.path.join(dir, pdb)
    cmd = "packmol < {} > {}".format(inpt_path, sf_path)
    os.system(cmd)
