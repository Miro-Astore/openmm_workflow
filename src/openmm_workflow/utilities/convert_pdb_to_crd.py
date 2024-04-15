import MDAnalysis as mda
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='topfile', help='Input topology file', required=True)
parser.add_argument('-ipdb', dest='pdbfile', help='Input coordinate pdb file',required=True)
parser.add_argument('-ocrd', dest='crdoutfile', help='Input coordinate pdb file', required=False)

args = parser.parse_args()

universe = mda.Universe(args.topfile,args.pdbfile)
if args.crdoutfile == None:
    out_file_name = args.pdbfile[:-4] + '.crd'
else:
    out_file_name = args.crdoutfile
    

writer_object = mda.coordinates.CRD.CRDWriter(out_file_name, extended=True)
writer_object.write(universe.atoms)
writer_object.close()
