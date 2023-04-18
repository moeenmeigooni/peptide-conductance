# 12 Apr 2022 Moeen Meigooni
#######################################

set sequence MPPPM ;# amino acid sequence
set sidelen 19.  ;# cube side length * 0.5
set n_conf 16

#######################################

for {set rep 0} {$rep < $n_conf} {incr rep} {

mol new ${sequence}_${rep}.psf
mol addfile ${sequence}_${rep}.pdb

package require solvate

# align S-S to Z axis
set all [atomselect top all]
set ss_pos [[atomselect top "name SD"] get {x y z}]
set ss_ax [vecnorm [vecsub [lindex $ss_pos 1] [lindex $ss_pos 0]]]
$all move [transvecinv $ss_ax]
$all move [transaxis y -90]
$all moveby [vecscale -1.0 [measure center $all]]
$all writepdb ${sequence}_align_${rep}.pdb

# solvate
mol new ${sequence}_${rep}.psf
mol addfile ${sequence}_align_${rep}.pdb
set all [atomselect top all]
set lims [measure minmax $all]
set maxs [vecsub "$sidelen $sidelen $sidelen" [lindex $lims 1]]
set mins [vecscale -1. [vecsub "-$sidelen -$sidelen -$sidelen" [lindex $lims 0]]]

solvate ${sequence}_${rep}.psf ${sequence}_align_${rep}.pdb -o ${sequence}_water_${rep} +x [lindex $maxs 0] +y [lindex $maxs 1] +z [lindex $maxs 2] -x [lindex $mins 0] -y [lindex $mins 1] -z [lindex $mins 2] 

mol new ${sequence}_water_${rep}.psf
mol addfile ${sequence}_water_${rep}.pdb

# translate coords for OpenMM-style PBC
[atomselect top all] moveby "$sidelen $sidelen $sidelen"
[atomselect top all] writepdb ${sequence}_water_${rep}.pdb

mol delete all

}

exit 

