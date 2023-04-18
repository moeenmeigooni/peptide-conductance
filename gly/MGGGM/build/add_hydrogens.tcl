# Moeen Meigooni 27 Apr 2022
#
#############################

set sequence MGGGM
set n_conf 16
set toppar_dir /Scr/meigoon2/toppar_c36_jul21

#############################
package require psfgen
resetpsf

##########
# make this operate on multiple conformers for {set i 0} {$i < $n} {incr i} etc

for {set rep 0} {$rep < $n_conf} {incr rep} {
  resetpsf
  topology $toppar_dir/top_all36_prot.rtf
  segment PEP { first NTER; last CTER; pdb ${sequence}_noh_${rep}.pdb }
  coordpdb ${sequence}_noh_${rep}.pdb PEP

  guesscoord

  writepsf ${sequence}_${rep}.psf
  writepdb ${sequence}_${rep}.pdb
}
exit

