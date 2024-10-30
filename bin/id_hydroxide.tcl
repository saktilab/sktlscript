proc do_highlight2 {args} {
   global highlight2 molid selid2
 #  set molid [molinfo top]
 #  set selid 2
   set frame [molinfo $molid get frame]
   if {[info exists highlight2($frame)]} then {
       mol modselect $selid2 $molid "$highlight2($frame)"
       set hom [atomselect top $highlight2($frame)]
       mol selupdate $hom 0 1
   }
}


set molid [molinfo top]
set selid2 3
set n [molinfo $molid get numframes]
set highlight2(0) none
set fp [open "id_hydroxide.dat" r]
for {set i 0} {$i < $n} {incr i} {
  set highlight2($i) [gets $fp]
#trajectory_path $h3op yellow 3 1
}
close $fp

trace variable vmd_frame($molid) w do_highlight2
animate goto start
