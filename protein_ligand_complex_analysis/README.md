## Protein–Ligand Interaction Analysis Workflow

This repository contains in-house scripts used to identify multiple classes of protein–ligand interactions from protein–ligand complex PDB files and to compile the detected interactions into a final interaction matrix.

## Start here

1. Place your input PDB complex files inside `Protein_ligand_complex_file/`
2. Run all commands from the project root directory which has the 'Protein_ligand_complex_file' folder
3. all Python scripts used are inside `scripts/`
4. See `Result` for a sample completed run

## Software requirements:
The scripts were executed in two Python environments. AromaticRings.py was run using /opt/miniconda2/bin/python (Python 2.7.18, Anaconda) with pybel and openbabel installed. The remaining analysis and matrix-processing scripts were run in a Python 3.8.10 environment; key dependencies included NumPy 1.20.2 and Pandas 1.2.4.

1. AromaticRings.py
   - Python 2.7.18 (Anaconda)
   - pybel
   - openbabel

2. Interaction scan scripts
   - Python 3.8.10
   - numpy 1.20.2

3. Matrix/filtering scripts
   - Python 3.8.10
   - numpy 1.20.2
   - pandas 1.2.4

## What this workflow does:

Given a folder of protein–ligand complex **PDB** files, this pipeline detects:

- binding-site contacts
- hydrogen bonds (ligand donor / ligand acceptor)
- chalcogen bonds
- halogen bonds
- aromatic interactions
- hydrophobic interactions

The resulting interaction files are then organized, cleaned, merged where needed, and converted into a final interaction matrix.

-------------------------------------------------------------------------------------------------------------------------------------------------------------

## The following commands are to be run :

1. Create the working directory with the 'Protein_ligand_complex_file' folder which contains the complex file/files

From the project directory:

mkdir Interactions
ls -1 Protein_ligand_complex_file/*.pdb > filelist.txt

2. Run the primary interaction-detection scripts

python scripts/Binding_site_prediction.py filelist.txt Interactions/Binding_sites.csv Interactions/Binding_sites_error.txt > Interactions/BindingSite.log 2>&1

python scripts/A.H-D_lig_donor_scan.py filelist.txt Interactions/Hbonds_donor.csv Interactions/Hbonds_donor_error.txt > Interactions/Hbonds_donor.log 2>&1

python scripts/A.H-D_lig_acceptor_scan.py filelist.txt Interactions/Hbonds_acceptor.csv Interactions/Hbonds_acceptor_error.txt > Interactions/Hbonds_acceptor.log 2>&1

python scripts/chalcogen_bond_lig_donor_scan.py filelist.txt Interactions/Chalcogen_donor.csv Interactions/Chalcogen_donor_error.txt > Interactions/Chalcogen_donor.log 2>&1

python scripts/chalcogen_bond_lig_acceptor_scan.py filelist.txt Interactions/Chalcogen_acceptor.csv Interactions/Chalcogen_acceptor_error.txt > Interactions/Chalcogen_acceptor.log 2>&1

python scripts/halogen_bond_scan.py filelist.txt Interactions/Halogen_interaction.csv Interactions/Halogen_interaction_error.txt > Interactions/Halogen_interaction.log 2>&1

These scripts generate:

output interaction tables
error files listing any PDBs that failed
log files

The binding-site script writes its output and error file from command-line arguments, and the hydrogen-bond donor script follows the same pattern.

3. Generate aromatic ring definitions for each complex

Create the ring-output directory:

mkdir AromaticRings

Then run:

for i in Protein_ligand_complex_file/*.pdb; do python scripts/AromaticRings.py "$i" done

This step identifies aromatic rings in each ligand and writes a corresponding .aromaticrings.txt file into AromaticRings/. AromaticRings.py uses Open Babel/Pybel and records ring index, ring type, ring size, aromaticity, and ring-atom coordinates.

4. Run aromatic interaction scans

Run all aromatic interaction scripts:

python scripts/Aromatic_interaction_f2f_5mem_scan.py > aromamtic_5mem_f2f.log 2>&1
python scripts/Aromatic_interaction_f2f_6mem_scan.py > aromamtic_6mem_f2f.log 2>&1
python scripts/Aromatic_interaction_f2f_7mem_scan.py > aromamtic_7mem_f2f.log 2>&1
python scripts/Aromatic_interaction_f2h_5mem_scan.py > aromamtic_5mem_f2h.log 2>&1
python scripts/Aromatic_interaction_f2h_6mem_scan.py > aromamtic_6mem_f2h.log 2>&1
python scripts/Aromatic_interaction_f2h_7mem_scan.py > aromamtic_7mem_f2h.log 2>&1

These scripts use the AromaticRings/*.aromaticrings.txt files generated in the previous step. The aromatic scripts explicitly open files from the AromaticRings/ directory and match ring atoms by ring size.

5. Run hydrophobic interaction scans

python scripts/hydrophobic_C_in_ligand_scan.py > hydrophobic_C_in_ligand.log 2>&1
python scripts/hydrophobic_C_in_protein_scan.py > hydrophobic_C_in_protein.log 2>&1

These generate raw hydrophobic interaction outputs that will be filtered later. This run order matches the command log from the original workflow.

6. Organize outputs into clean subfolders

Move into the Interactions/ folder and create subfolders:

cd Interactions

mkdir ERROR
mv ../error* ERROR/ 2>/dev/null || true
mv Binding_sites_error.txt Hbonds_acceptor_error.txt Hbonds_donor_error.txt ERROR/ 2>/dev/null || true

mkdir LOG
mv *.log LOG/ 2>/dev/null || true
mv ../*.log LOG/ 2>/dev/null || true

mkdir AromaticInteractions
mv ../out_aromatic_interaction_f2* AromaticInteractions/

mkdir Hydrophobic
mv ../out_hydrophobic_C_in_* Hydrophobic/

mkdir FinalInteractionFiles
mv *.csv FinalInteractionFiles/
mv ../out_chalcogen_bond_lig_* FinalInteractionFiles/
mv ../out_halogen_bond_scan.txt FinalInteractionFiles/

This step is only for cleaning and organizing the outputs. It separates:

error files
logs
aromatic interaction files
hydrophobic interaction files
final curated interaction files

7. Deduplicate and combine aromatic interaction outputs
Move into the aromatic interaction folder:
cd AromaticInteractions

Keep the first 9 tab-separated fields, sort the records, and remove duplicates:

cut -f1-9 out_aromatic_interaction_f2f_5mem_scan.txt | sort | uniq > out_aromatic_interaction_f2f_5mem_scan_uniq.txt
cut -f1-9 out_aromatic_interaction_f2f_6mem_scan.txt | sort | uniq > out_aromatic_interaction_f2f_6mem_scan_uniq.txt
cut -f1-9 out_aromatic_interaction_f2f_7mem_scan.txt | sort | uniq > out_aromatic_interaction_f2f_7mem_scan_uniq.txt
cut -f1-9 out_aromatic_interaction_f2h_5mem_scan.txt | sort | uniq > out_aromatic_interaction_f2h_5mem_scan_uniq.txt
cut -f1-9 out_aromatic_interaction_f2h_6mem_scan.txt | sort | uniq > out_aromatic_interaction_f2h_6mem_scan_uniq.txt
cut -f1-9 out_aromatic_interaction_f2h_7mem_scan.txt | sort | uniq > out_aromatic_interaction_f2h_7mem_scan_uniq.txt

Combine ring-size-specific outputs into one non-redundant file for each aromatic interaction class:

cat out_aromatic_interaction_f2f_5mem_scan_uniq.txt \
    out_aromatic_interaction_f2f_6mem_scan_uniq.txt \
    out_aromatic_interaction_f2f_7mem_scan_uniq.txt | sort | uniq > aromatic_interaction_f2f_combined.txt

cat out_aromatic_interaction_f2h_5mem_scan_uniq.txt \
    out_aromatic_interaction_f2h_6mem_scan_uniq.txt \
    out_aromatic_interaction_f2h_7mem_scan_uniq.txt | sort | uniq > aromatic_interaction_f2h_combined.txt
    
Copy the combined aromatic files to the final folder:

cp aromatic_interaction_f2* ../FinalInteractionFiles/
cd ..

In this step:

cut -f1-9 keeps only the first 9 tab-separated columns
sort | uniq removes duplicates
cat merges the 5-, 6-, and 7-member ring outputs into one final combined file for each aromatic interaction type

This is exactly the post-processing performed in the original terminal workflow.

8. Filter hydrophobic interaction outputs

Move into the hydrophobic directory:

cd Hydrophobic

python ../scripts/filter.py

Copy the filtered outputs into the final folder:

cp *_scan_filtered.txt ../FinalInteractionFiles/
cd ..

9. Build the final interaction matrix
Inside Interactions/, copy the file list as A:

cp ../filelist.txt A

Then Build the final interaction matrix:

Before running this step, edit `scripts/make_matrix_bindingsite.py` and `scripts/filter_matrix_count.py` to include the binding-pocket residues of your target protein.

In the current version, these scripts are configured for the DPP4 (PDB: 1X70) binding-pocket residues:
- Glu205
- Glu206
- Ser209
- Phe357
- Arg358
- Tyr547
- Tyr662

Also update the matrix/output filenames if needed (for example, `matrix_DPP4_Chitraline.csv`). The current scripts are hardcoded for these residue definitions and filenames. 

After editing, copy the scripts from the `scripts/` folder into the `Interactions/` directory and run:

cp ../scripts/make_matrix_bindingsite.py .
cp ../scripts/filter_matrix_count.py .

python make_matrix_bindingsite.py > matrix_bindingsite.log 2>&1
python filter_matrix_count.py

This step generates three tab-separated output files:

make_matrix_DPP4_Chitraline.csv
Intermediate file containing, for each ligand and each interaction file, whether the ligand interacts with each of the seven binding-pocket residues.
matrix_DPP4.csv
Final wide interaction matrix containing one row per ligand and residue-wise boolean columns for each interaction class.
matrix_DPP4_interactions_count.csv
Summary file that collapses all interaction classes and reports whether each of the seven residues is contacted by any interaction, along with AllInteractions, the total number of binding-pocket residues contacted by that ligand.
