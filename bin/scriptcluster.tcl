set sel [atomselect top all]
set cenxyz [measure center $sel]
set cenx [lindex $cenxyz 0]
set ceny [lindex $cenxyz 1]
set cenz [lindex $cenxyz 2]

for {set i 1} {$i<=20} {incr i} {
 set thres [expr 0.5*$i]
 set sel [atomselect top "(x-$cenx)**2 + (y-$ceny)**2 + (z-$cenz)**2 < $thres**2"]
 puts [$sel num]
 if {[$sel num]>0} {break}
}
set cenfrag [lindex [$sel get fragment] 0]
puts "same fragment as {within 3.5 of fragment $cenfrag}"
