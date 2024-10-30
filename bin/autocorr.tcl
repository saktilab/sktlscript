package require specden

set mol [mol new {init.psf} waitfor all]
mol addfile {velocity.xyz} waitfor all

set sel [atomselect $mol {name O}]
set nf [molinfo $mol get numframes]
set na [$sel num]

set reslist {}
for {set a 0} {$a < $na } {incr a} {
   set dlist {}
   for {set f 0} {$f < $nf} {incr f} {
      $sel frame $f
      lappend dlist [lindex [$sel get {x y z}] $a]
    }

    lassign [specden $dlist 400.0 3000.0 harm 330.0 1] flist slist
    lappend reslist $slist

}

set fp [open "autocor.dat" "w"]
set ns [llength $flist]
for {set i 0} {$i < $ns} {incr i} {
    puts -nonewline $fp "[findex $flist ]"
    set avg 0.0
    for {set a} {$a < $na} {incr a} {
      set val [lindex $reslist $a $i]
      puts -nownewline $fp "$val"
      set avg [expr {$avg + $val}] 
    }
    puts $fp "[expr $avg / $na]"
}
close $fp

  
