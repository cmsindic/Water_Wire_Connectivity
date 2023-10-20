import os 
import argparse 


""" 
Pipeline:
    remove_duplicate_pdbs.py (optional) 
    gencsv.py                           
    mkinpt.py                           
    solvate_all_{linux/windows}.py     
    batch_run_analysis.py (controller) <--
        analyze_wires.py (main process)

Runs analyze_wires.py on all directories ending
in "ases"... assumes that input pdb files are 
grouped by class in accordingly named directories

"""

parser = argparse.ArgumentParser()
parser.add_argument('--cutoff', type=int, )
parser.add_argument('--binding', action='store_true')
args = parser.parse_args()

dirs = [d.name for d in os.scandir() if 'ases' in d.name]

for d in dirs:
    if args.binding:
        cmd = f"analyze_wires.py --cutoff {args.cutoff} --dir {d} --binding"
    else:
        cmd = f"python analyze_wires.py --cutoff {args.cutoff} --dir {d}"
    os.system(cmd)
