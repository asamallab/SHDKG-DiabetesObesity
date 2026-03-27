prot_aids = {}
lig_aids = {}

for i in open('../FinalInteractionFiles/Hbonds_acceptor.csv'):
    tmp = [j.strip() for j in i.strip().split('\t')]
    fil = tmp[0].strip()
    laid = [tmp[4]]
    paid = [tmp[11]]
    try:
        prot_aids[fil] += paid
        lig_aids[fil] += laid
    except KeyError:
        prot_aids[fil] = paid
        lig_aids[fil] = laid
for i in open('../FinalInteractionFiles/Hbonds_donor.csv'):
    tmp = [j.strip() for j in i.strip().split('\t')]
    fil = tmp[0].strip()
    laid = [tmp[11]]
    paid = [tmp[4]]
    try:
        prot_aids[fil] += paid
        lig_aids[fil] += laid
    except KeyError:
        prot_aids[fil] = paid
        lig_aids[fil] = laid
for i in open('../FinalInteractionFiles/out_chalcogen_bond_lig_acceptor_scan.txt'):
    tmp = [j.strip() for j in i.strip().split('\t')]
    fil = tmp[0].strip()
    laid = [tmp[4]]
    paid = [tmp[10],tmp[12]]
    try:
        prot_aids[fil] += paid
        lig_aids[fil] += laid
    except KeyError:
        prot_aids[fil] = paid
        lig_aids[fil] = laid
for i in open('../FinalInteractionFiles/out_chalcogen_bond_lig_donor_scan.txt'):
    tmp = [j.strip() for j in i.strip().split('\t')]
    fil = tmp[0].strip()
    laid = [tmp[10],tmp[12]]
    paid = [tmp[4]]
    try:
        prot_aids[fil] += paid
        lig_aids[fil] += laid
    except KeyError:
        prot_aids[fil] = paid
        lig_aids[fil] = laid
for i in open('../FinalInteractionFiles/out_halogen_bond_scan.txt'):
    tmp = [j.strip() for j in i.strip().split('\t')]
    fil = tmp[0].strip()
    laid = [tmp[11],tmp[14]]
    paid = [tmp[4], tmp[7]]
    try:
        prot_aids[fil] += paid
        lig_aids[fil] += laid
    except KeyError:
        prot_aids[fil] = paid
        lig_aids[fil] = laid

print (len(prot_aids), len(lig_aids))
#print (prot_aids['NSP15_DOCKED_LIGAND_COMBINED/NSP15_100005_ligand_1_combined.pdb'])
#print (lig_aids['NSP15_DOCKED_LIGAND_COMBINED/NSP15_100005_ligand_1_combined.pdb'])

fout1 = open('out_hydrophobic_C_in_ligand_scan_filtered.txt','w')
for i in open('out_hydrophobic_C_in_ligand_scan.txt'):
    tmp = i.strip().split('\t')
    fil = tmp[0].strip()
    laid = tmp[9].strip()
    paid = tmp[4].strip()
    if fil in prot_aids.keys():
        if laid in lig_aids[fil] or paid in prot_aids[fil]:
            continue
        else:
            fout1.write(i)
    else:
        fout1.write(i)
fout2 = open('out_hydrophobic_C_in_protein_scan_filtered.txt','w')
for i in open('out_hydrophobic_C_in_protein_scan.txt'):
    tmp = i.strip().split('\t')
    fil = tmp[0].strip()
    laid = tmp[9].strip()
    paid = tmp[4].strip()
    if fil in prot_aids.keys():
        if laid in lig_aids[fil] or paid in prot_aids[fil]:
            continue
        else:
            fout2.write(i)
    else:
        fout2.write(i)
fout1.close()
fout2.close()

