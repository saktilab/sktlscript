proc make_rotation_animated_gif {} {
	set frame 0
	for {set i 0} {$i < 360} {incr i 5} {
		set filename snap.[format "%04d" $frame].rgb
		render snapshot $filename
		incr frame
		rotate y by 20
	}
	exec convert -delay 30 -loop 4 snap.*.rgb movie.gif
	rm snap.*.rgb
}
