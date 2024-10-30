# Overlays
 # Author: Jordi Cohen
 #
 # For now, implements a grid overlay
 namespace eval ::Overlays:: {
   set grid_on 0
   set grid_dirty 1  ;# ruler needs a refresh
   set grid_scale 0
   set grid_scale_graphics_id -1
   set grid_mol -1
 }
 proc ::Overlays::setup_grid {} {
   variable grid_mol
   variable grid_on
   variable grid_dirty
   set top [molinfo top]
   set grid_mol [mol new]
   mol rename $grid_mol "Grid"
   if {$top >= 0} {
     mol top $top
     molinfo $grid_mol set scale_matrix [molinfo $top get scale_matrix]
   }
   redraw_grid
   trace add variable ::vmd_logfile write ::Overlays::logfile_cb
   set grid_on 1
   set grid_dirty 1
 }
 proc ::Overlays::remove_grid {} {
   variable grid_mol
   variable grid_on
   trace remove variable ::vmd_logfile write ::Overlays::logfile_cb
   mol delete $grid_mol
   set grid_on 0
 }
 proc ::Overlays::redraw_grid {} {
   variable grid_mol
   variable grid_on
   variable grid_scale
   variable grid_dirty
   variable grid_scale_graphics_id
   molinfo $grid_mol set center_matrix [list [transidentity]]
   molinfo $grid_mol set rotate_matrix [list [transidentity]]
   molinfo $grid_mol set global_matrix [list [transidentity]]
   set realscale [lindex [molinfo $grid_mol get scale_matrix] 0 0 0]
   set scale [expr round(-log10($realscale))-1]
   if {$scale != $grid_scale} {set grid_dirty 1}
   set display_ratio [expr 1.*[lindex [display get size] 0]/[lindex [display get
 size] 1]]
   set div [expr pow(10,$scale)]
   if {$grid_dirty} {
     set grid_scale $scale
     set minx [expr -200*$div]
     set maxx [expr 200*$div]
     set miny [expr -50*$div]
     set maxy [expr 50*$div]
     graphics $grid_mol delete all
     graphics $grid_mol  color gray
   #  draw material Transparent
     for {set tick $minx} {$tick <= $maxx} {set tick [expr $tick + $div]} {
       graphics $grid_mol  line [list $tick $miny 0] [list $tick $maxy 0] width 1
 style dashed
     }
     for {set tick $miny} {$tick <= $maxy} {set tick [expr $tick + $div]} {
       graphics $grid_mol  line [list $minx $tick 0] [list $maxx $tick 0] width 1
 style dashed
     }
     draw color gray
 #  draw material Opaque
     for {set tick $minx} {$tick <= $maxx} {set tick [expr $tick + 10.*$div]}
 {
       graphics $grid_mol line [list [expr $tick] $miny 0] [list [expr $tick]
 $maxy 0] width 2
     }
     for {set tick $miny} {$tick <= $maxy} {set tick [expr $tick + 10.*$div]}
 {
       graphics $grid_mol  line [list $minx $tick 0] [list $maxx $tick 0] width 2
     }
     set grid_scale_graphics_id [graphics $grid_mol text [list [expr
 1.2*$display_ratio/$realscale] [expr -1.4/$realscale] 0] "[format
 "%g" $div]A" size 0.8]
     set grid_dirty 0
   } else {
     graphics $grid_mol delete $grid_scale_graphics_id
     set grid_scale_graphics_id [graphics $grid_mol text [list [expr
 1.2*$display_ratio/$realscale] [expr -1.4/$realscale] 0] "[format
 "%g" $div]A" size 0.8]
   }
 }
 proc ::Overlays::logfile_cb { args } {
   variable grid_mol
   # Check for display transforms
   if { [string match "rotate*" $::vmd_logfile] || [string match
 "translate*" $::vmd_logfile] || \
             [string match "scale*" $::vmd_logfile] || [string match
 "display*" $::vmd_logfile]} {
     redraw_grid
   }
 }
 proc ::Overlays::overlay {args} {
   variable grid_mol
   variable grid_on
   set overlay [lindex $args 0]
   set state [string is true [lindex $args 1]]
   if {"$overlay" == "grid"} {
     if {$state} {
       if {$grid_on} {remove_grid}
       setup_grid
     } else {
       if {$grid_on} {remove_grid}
     }
   } else {
     puts "overlay: Unknown overlay"
   }
 }
 proc overlay {args} {
   eval ::Overlays::overlay $args
 }
 overlay grid on
