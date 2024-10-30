proc do_highlight_max_1 {args} {
   global highlight_max_1 molid selid_max_1
 #  set molid [molinfo top]
 #  set selid 2
   set frame [molinfo $molid get frame]
   if {[info exists highlight_max_1($frame)]} then {
       mol modselect $selid_max_1 $molid "$highlight_max_1($frame)"
       set h3op [atomselect top $highlight_max_1($frame)]
       mol selupdate $h3op 0 1
   }
}


set molid [molinfo top]
set selid_max_1 [molinfo top get numreps]
mol addrep 0
mol modstyle $selid_max_1 0 VDW 1.0 12
mol modcolor $selid_max_1 0 ColorId 4
mol modmaterial $selid_max_1 0 Transparent
set n [molinfo $molid get numframes]
set highlight_max_1(0) none
set fp [open "id_hydronium1.dat" r]
for {set i 0} {$i < $n} {incr i} {
  set highlight_max_1($i) [gets $fp]
}
close $fp

trace variable vmd_frame($molid) w do_highlight_max_1
animate goto start
