set num [molinfo top get numframes]
set ox [atomselect top {name O}]
set all [atomselect top {all}]

foreach i [$ox get index] {
    set sel($i) [atomselect top "exwithin 1.30 of index $i"]
}


for {set n 0} {$n < $num} {incr n} {
   set bc {}
   foreach i [$ox get index] {
      $sel($i) frame $n
      $sel($i) update
      $all frame $n
      $all set user 0
      lappend bc [$sel($i) num]
   }
   $ox frame $n
   $ox set user $bc
   unset bc
}

foreach i [$ox get index] {
   $sel($i) delete
}

$ox delete
$all delete
unset ox all sel i n
