#!/usr/bin/env python3
# Updated for 1X70 binding-pocket residues:
# Glu205, Glu206, Ser209, Phe357, Arg358, Tyr547, Tyr662
# Output matrix file: matrix_1X70.csv

import os
import pandas as pd

# Ligand IDs (one per line) in file "A"
files = [i.strip() for i in open('A') if i.strip()]

outfiles1 = [
    'FinalInteractionFiles/Binding_sites.csv',
    'FinalInteractionFiles/out_halogen_bond_scan.txt',
    'FinalInteractionFiles/Hbonds_donor.csv',
    'FinalInteractionFiles/out_hydrophobic_C_in_protein_scan_filtered.txt',
    'FinalInteractionFiles/out_hydrophobic_C_in_ligand_scan_filtered.txt',
    'FinalInteractionFiles/aromatic_interaction_f2h_combined.txt',
    'FinalInteractionFiles/out_chalcogen_bond_lig_donor_scan.txt'
]

outfiles2 = [
    'FinalInteractionFiles/Hbonds_acceptor.csv'
]

def ensure_nonempty(path: str, ncols: int = 9) -> None:
    """
    If an interaction file is empty (0 bytes/0 lines), pandas.read_csv() crashes.
    This makes a harmless placeholder row with enough tab columns so indexing works.
    """
    try:
        if (not os.path.exists(path)) or os.path.getsize(path) == 0:
            row = ["PLACEHOLDER"] + ["NA"] * (ncols - 1)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                f.write("\t".join(row) + "\n")
    except OSError:
        # If something odd happens with filesystem, let pandas raise a clear error later.
        pass

# ---- STEP 1: Build per-file, per-ligand residue-hit booleans ----
with open('make_matrix_DPP4_Chitraline.csv', 'w') as fout:
    # outfiles1: residue is assumed in columns [2]=RESNAME, [3]=RESNUM
    for i in outfiles1:
        print(i)
        ensure_nonempty(i, ncols=9)
        df = pd.read_csv(i, header=None, sep='\t', dtype=str)

        for j in files:
            tmp1 = df[(df[0] == j)]

            chk1 = tmp1[(tmp1[2] == 'GLU') & (tmp1[3] == '205')]
            chk2 = tmp1[(tmp1[2] == 'GLU') & (tmp1[3] == '206')]
            chk3 = tmp1[(tmp1[2] == 'SER') & (tmp1[3] == '209')]
            chk4 = tmp1[(tmp1[2] == 'PHE') & (tmp1[3] == '357')]
            chk5 = tmp1[(tmp1[2] == 'ARG') & (tmp1[3] == '358')]
            chk6 = tmp1[(tmp1[2] == 'TYR') & (tmp1[3] == '547')]
            chk7 = tmp1[(tmp1[2] == 'TYR') & (tmp1[3] == '662')]

            fout.write('\t'.join([
                i, j,
                str(len(chk1) != 0), str(len(chk2) != 0), str(len(chk3) != 0),
                str(len(chk4) != 0), str(len(chk5) != 0), str(len(chk6) != 0),
                str(len(chk7) != 0)
            ]) + '\n')

    # outfiles2: residue is assumed in columns [7]=RESNAME, [8]=RESNUM
    for k in outfiles2:
        print(k)
        ensure_nonempty(k, ncols=9)
        df2 = pd.read_csv(k, header=None, sep='\t', dtype=str)

        for l in files:
            tmp2 = df2[(df2[0] == l)]

            chk1_2 = tmp2[(tmp2[7] == 'GLU') & (tmp2[8] == '205')]
            chk2_2 = tmp2[(tmp2[7] == 'GLU') & (tmp2[8] == '206')]
            chk3_2 = tmp2[(tmp2[7] == 'SER') & (tmp2[8] == '209')]
            chk4_2 = tmp2[(tmp2[7] == 'PHE') & (tmp2[8] == '357')]
            chk5_2 = tmp2[(tmp2[7] == 'ARG') & (tmp2[8] == '358')]
            chk6_2 = tmp2[(tmp2[7] == 'TYR') & (tmp2[8] == '547')]
            chk7_2 = tmp2[(tmp2[7] == 'TYR') & (tmp2[8] == '662')]

            fout.write('\t'.join([
                k, l,
                str(len(chk1_2) != 0), str(len(chk2_2) != 0), str(len(chk3_2) != 0),
                str(len(chk4_2) != 0), str(len(chk5_2) != 0), str(len(chk6_2) != 0),
                str(len(chk7_2) != 0)
            ]) + '\n')

outfiles = outfiles1 + outfiles2

# ---- STEP 2: Pivot into one wide matrix per ligand ----
with open('matrix_DPP4.csv', 'w') as fout2:
    header = ['Filename']
    for i in [
        'BindingSite', 'Halogen_bond', 'Hbond_donor', 'Hydrophobic_CinProtein',
        'Hydrophobic_CinLigand', 'Aromatic_f2h', 'chalcogen_doner', 'Hbond_acceptor'
    ]:
        for j in ['Glu205', 'Glu206', 'Ser209', 'Phe357', 'Arg358', 'Tyr547', 'Tyr662']:
            header.append(j + '_' + i)
    fout2.write('\t'.join(header) + '\n')

    # Read step1 once into memory (still simple, but faster than reopening for every k)
    step1_lines = [line.strip().split('\t') for line in open('make_matrix_DPP4_Chitraline.csv') if line.strip()]

    for j in files:
        mat = [j]
        for k in outfiles:
            found = False
            for tmp in step1_lines:
                outf = tmp[0]
                ff = tmp[1]
                chk = tmp[2:]
                if len(chk) != 7:
                    print('Error', j, k, tmp)
                    break
                if ff == j and outf == k:
                    mat += chk
                    found = True
                    break
            if not found:
                # If an interaction file has no entries for this ligand, fill with 7 False
                mat += ['False'] * 7

        fout2.write('\t'.join(mat) + '\n')

