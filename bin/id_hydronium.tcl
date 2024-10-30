proc do_highlight {args} {
   global highlight molid selid
 #  set molid [molinfo top]
 #  set selid 2
   set frame [molinfo $molid get frame]
   if {[info exists highlight($frame)]} then {
       mol modselect $selid $molid "$highlight($frame)"
       set h3op [atomselect top $highlight($frame)]
       mol selupdate $h3op 0 1
   }
}


set molid [molinfo top]
set selid 2
set n [molinfo $molid get numframes]
set highlight(0) none
set fp [open "id_hydronium.dat" r]
for {set i 0} {$i < $n} {incr i} {
  set highlight($i) [gets $fp]
#trajectory_path $h3op yellow 3 1
}
close $fp

trace variable vmd_frame($molid) w do_highlight
animate goto start
