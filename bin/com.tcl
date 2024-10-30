set xmol [mol load psf gblmix.psf]
mol addfile traject.xyz waitfor all $xmol
pbc set {21.7291788807 21.7291788807 21.7388403991} -all
  # pbc wrap -all
# pbc join bonded -now
# pbc wrap -now
# pbc unwrap -first now

set ntracked 66
set ymol [mol new atoms $ntracked];
set mov ""
for {set i 0} {$i < $ntracked} {incr i} {lappend mov [atomselect $ymol "index $i"]}
set all [atomselect $ymol all]
$all set radius 1.0
$all delete
animate write psf com_traj.psf $ymol

set sol [atomselect $xmol all]
set resids [lsort -unique [$sol get resid]]
$sol delete
set sel ""
foreach resid $resids { lappend sel [atomselect $xmol "all and resid $resid"] }

for {set i 0} {$i < [molinfo $xmol get numframes]} {incr i} {
animate dup $ymol
set atomindex 0
foreach selection $sel {
  $selection frame $i
  set center [measure center $selection weight mass]
  [lindex $mov $atomindex] moveto $center
  incr atomindex 
}
}
animate write xyz com_traj.xyz $ymol

foreach selection [atomselect list] {$selection delete}