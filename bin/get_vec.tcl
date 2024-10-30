
set nf [molinfo top get numframes]
set outfile [open "dist.dat" w]

set all [atomselect top all]
set no3 [atomselect top "index 1500"]
for {set i 0} {$i < $nf} {incr i} {
$all frame $i
$no3 frame $i
set x [measure center $all]
set y [measure center $no3]
set length [vecdist $x $y]
$all update
$no3 update
   puts $outfile "$i $length"   

}

close $outfile
