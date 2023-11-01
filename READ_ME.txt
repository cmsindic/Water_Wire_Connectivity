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
    
