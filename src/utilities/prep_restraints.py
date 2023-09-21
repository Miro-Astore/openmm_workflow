import MDAnalysis as mda
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='topfile', help='Input topology file', required=True)
parser.add_argument('-ipdb', dest='pdbfile', help='Input coordinate pdb file',required=True)
args = parser.parse_args()

universe = mda.Universe(args.topfile, args.pdbfile)

backbone = universe.select_atoms('backbone and not name H*')
sidechains = universe.select_atoms('protein and (not backbone and not name H*)')

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
