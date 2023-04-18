#!/bin/bash
###########################
# user defined var
PEPTIDE_SEQUENCE='MPPPM'
###########################

# make peptide w/o hydrogens
/Scr/meigoon2/miniconda3/bin/python make_peptides.py $PEPTIDE_SEQUENCE
# add hydrogens
vmd -dispdev text -e add_hydrogens.tcl
# align and solvate
vmd -dispdev text -e solvate.tcl



