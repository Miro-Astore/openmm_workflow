import MDAnalysis as mda
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='topfile', help='Input topology file', required=True)
parser.add_argument('-ipdb', dest='pdbfile', help='Input coordinate pdb file',required=True)
parser.add_argument('-bb', dest='backbone_selection', help='Selection string used to find backbone atoms.', default='protein and backbone and not name H*')
parser.add_argument('-sc', dest='sidechain_selection', help='Selection string used to find sidechain atoms.', default='protein and not backbone and not name H*')
args = parser.parse_args()

universe = mda.Universe(args.topfile, args.pdbfile)

backbone = universe.select_atoms(args.backbone_selection)
sidechains = universe.select_atoms(args.sidechain_selection)

try:
        os.mkdir('restraints')
except FileExistsError:
        pass

with open ('restraints/prot_pos.txt','w') as out_file:
    for i in backbone.atoms:
        out_file.write(str(i.index) + ' BB')

    for i in sidechains.atoms:
        out_file.write(str(i.index) + ' SC')

out_file.close()
