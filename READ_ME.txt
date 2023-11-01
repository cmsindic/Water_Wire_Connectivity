This series of scripts is designed to take pdb files, located
in directories in accordance with their enzyme class, and 
solvate them using PACKMOL before assessing for whether water
can be found within a cutoff distance of their active site 
residues. The number of water molecules within this cutoff 
distance of the active site residues are then treated as a 
method of seeding water wires, water within cutoff distance 
of one another. Waters in each enzyme model are then recursively
added the the water wire emanating from its active site. 

Moreover, as this series of scripts is a method for testing 
the contiguity of active sites and bulk solvent in enzyme models,
for each model, the number of bulk solvent molecules, which are 
further than the set cutoff distance from the enzyme atoms, that 
are within this cutoff distance from water wire molecules are 
counted. 

PIPELINE:

gencsv.py                         

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


mkinpt.py 
        
    This code takes csv files made by gencsv.py
    and converts them into inp files for 
    PACKMOL to use to solvate the PDB files 
    corresponding with each csv file. 
    
    It requires model.inp to be in the working
    directory to be used as a template for its 
    output.


solvate_all_{linux/windows}.py

    Runs packmol on all inp/pdb file pairs in ./*ases
    that have not already been solvated. 


analyze_wires.py

    This program takes pdb files {ID}_solvated.pdb 
    that have been solvated by packmol, extracts active
    site locations from corresponding .index files which
    contain active site line numbers in the ORIGINAL pdb
    file {ID}.pdb, writes the active site lines into 
    their own pdb, then reads this pdb to obtain active 
    site residues in the solvated file. 
    
    Both the solvated file and the temporary active site
    pdb file are parsed with Bio.PDB.PDBParser. The 
    program then establishes water in the solvated pdb 
    as SOLVENT or NEAR ENZYME. It then forms a wire of 
    waters within distance CUTOFF of the active site atoms
    and propagates this wire, finding more and more waters
    within CUTOFF of each other. It then asks how many 
    solvent water molecules are within CUTOFF of any 
    component of the wire. This and more data are written 
    to a file in the target directory, 
        {ID}_wireInfo_cutoff_{CUTOFF}.csv.
    
