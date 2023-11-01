import os 

dirs = [d for d in os.scandir() if "ases" in d.name]

def get_model_id(file):
    return os.path.split(file.path)[-1][:-4]

     
for d in dirs:
    files = [f for f in os.scandir(d)]
    s = lambda x: "solvated" in x.name
    p = lambda x: "pdb" in x.name
    
    models = []
    for f in files:
        if len(f.name) != 8:
            continue
        if not p(f):
            continue
        if s(f):
            continue
        m = get_model_id(f)
        models.append(m)
    
    original_pdbs = [m + ".pdb" for m in models]
    for f in os.scandir(d):
        if not f.name in original_pdbs:
            os.remove(f)
        
                
            
        
        