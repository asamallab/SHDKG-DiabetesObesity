# Sat May  9 10:45:04 IST 2020
# vivekananth@imsc.res.in
# Program to identify the Aromatic Rings.
# input a pdb file of the ligand (First use vina_split to take the first model and use mgl scripts pdbqt_to_pdb.py to comvert it to pdb format,then use it as input here)

import sys
import pybel
import openbabel

fout = open('AromaticRings/'+sys.argv[1].strip().split('/')[-1].replace('.pdb','.aromaticrings.txt'),'w')

for mol in pybel.readfile( "pdb", sys.argv[1]):
    rings = str(len(mol.OBMol.GetSSSR()))
    ar=0
    for ind,ring in enumerate(mol.OBMol.GetSSSR()):
        rsize = str(ring.Size())
        aromatic = ring.IsAromatic()
        rtype = ring.GetType()
        if aromatic == True:
            ar=ar+1
            for a in mol.atoms:
                aidx = a.idx
                if ring.IsMember(a.OBAtom):
                    ax = str(a.OBAtom.GetX())
                    ay = str(a.OBAtom.GetY())
                    az = str(a.OBAtom.GetZ())
                    fout.write('\t'.join([str(ind+1),rtype,rsize,str(aromatic),ax,ay,az]) + '\n')
                else:
                    continue
        else:
            continue

fout.write('# Number of rings in the molecule: {}\n'.format(rings))
fout.write('# Number of Aromatic rings in the molecule: {}\n'.format(str(ar)))

fout.close()
