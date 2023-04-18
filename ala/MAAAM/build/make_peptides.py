# Moeen Meigooni 27 Apr 2022
# creates 
# usage: 
# >> python make_peptide.py KLVFF
# arg var can be any AA sequence

import PeptideBuilder as pb
import PeptideBuilder.Geometry as geo
import Bio.PDB
import numpy as np
from sys import argv

script, sequence = argv

def random_angle():
    # returns random value between -180, 180
    return (np.random.random() * 360.) - 180.

# define number of conformations to generate
n_conf = 16

for confid in range(n_conf):
    # use random phi, psi angles for each amino acid in the peptide
    geos = []
    for aa in sequence:
        g = geo.geometry(aa)
        g.phi = random_angle()
        g.psi_im1 = random_angle()
        geos.append(g)
    
    # build 
    structure = pb.initialize_res(geos[0])
    for i in range(1,len(sequence)):
        pb.add_residue(structure, geos[i])
    pb.add_terminal_OXT(structure)
    
    # save
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save('%s_noh_%i.pdb'%(sequence, confid))






