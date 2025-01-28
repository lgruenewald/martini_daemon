# Advanced visualization script for martini_daemon

# decodes a int32 list that was written to a file, returns a Tcl list of ints
proc decode_list {filename} {
	set fd [open $filename r]
	fconfigure $fd -translation binary
	set binary [read $fd]
	close $fd

	set result {}
	set length [string length $binary]

	for {set i 0} {$i < $length} {incr i 4} {
		set chunk [string range $binary $i [expr {$i + 3}]]
		binary scan $chunk i value
		lappend result $value
	}

	return $result
}

# decodes a list of lists, of variable length that was written to a file
# i32 of -1 (FFFFFFFF in hex) is used as a list separator, so this is not a
# valid value the lists can contain
proc decode_nested_list {filename} {
	set fd [open $filename r]
	fconfigure $fd -translation binary
	set binary [read $fd]
	close $fd

	set result {}
	set length [string length $binary]

	set i 0
	set c {}

	while {$i < $length} {
		set chunk [string range $binary $i [expr {$i + 3}]]
		binary scan $chunk i value

		if {$value == -1} {
			lappend result $c
			set c {}
		} else {
			lappend c $value
		}
		incr i 4
	}
	lappend result $c

	return $result
	
}

# opens a .gro and .xtc with CPK that can be used for visualizing e.g. the bonds
# removes the extra frame from the .gro
# supply molid if opening as a different molid, though the supplied molid
# should be manually verified to be the next one to be added
proc daemon_open {gro xtc {molid 0}} {
	mol new $gro
	mol addfile $xtc
	animate delete  beg 0 end 0 skip 0 $molid
	mol modstyle 0 $molid CPK 3.0 1.5 12.0 12.0
	animate pause
}

# colors molecule molid according to the daemon reporter output file npy
proc daemon_color {npy {molid 0}} {
	
	set dsel [ atomselect $molid all ]
	set all_frames [ decode_list $npy ]
	set len [ llength $all_frames ]
	# number of atoms, assumed to be costant
	set n [ $dsel num ]
	set n_frames [expr $len / $n]
	
	for {set frame 0} {$frame < $n_frames} {incr frame 1} {
		# for every frame...
		$dsel frame $frame
		set frame_start [expr $n * $frame]
		set frame_end [expr $frame_start + $n - 1]
		# lrange is inclusive on both ends
		set frame_data [ lrange $all_frames $frame_start $frame_end]
		# set custom data
		$dsel set user $frame_data
	}

	# if called multiple times after itself it should always clean up the 
	# reprsentations of the previous
	for {set i 0} {$i < 32} {incr i 1} {
		mol delrep 0 $molid
	}

	# this is how discrete colors are achieved
	for {set i 0} {$i < 32} {incr i 1} {
		mol color ColorID $i
		mol representation CPK 3.0 1.5 12.0 12.0
		mol selection user % 32 == $i
		mol material Opaque
		mol addrep $molid
	}

	# display some diagnostic data
	puts "daemon_color: loaded $len integers; atom n $n; frame n $n_frames;"
	
}

# displays bonds in frame $frame for molecule $molid according to the
# daemon bonds reporter output file $npy
proc daemon_bonds {npy frame {molid 0}} {
	set dsel [ atomselect $molid all ]
	set all_frames [ decode_nested_list $npy ]
	set len [ llength $all_frames ]
	set n [ $dsel num ]
	set n_frames [ expr $len / $n ]

	# sadly setbonds seems to not be frame specific (unlike User)
	$dsel frame $frame
	set frame_start [expr $n * $frame]
	set frame_end [expr $frame_start + $n - 1]
	set frame_data [lrange $all_frames $frame_start $frame_end]
	$dsel setbonds $frame_data
	animate goto $frame

	puts "daemon_bonds: loaded $len lists; atom n $n; frame n $n_frames;"
	puts "selected frame $frame data $frame_start to $frame_end"
}

proc render_bonds {npy {molid 0}} {
	set dsel [ atomselect $molid all ]
	set all_frames [ decode_nested_list $npy ]
	set len [ llength $all_frames ]
	set n [ $dsel num ]
	set n_frames [ expr $len / $n ]
	
	puts "render_bonds: loaded $len lists; atom n $n; frame n $n_frames;"
	for {set frame 0} {$frame < $n_frames} {incr frame 1} {
		$dsel frame $frame
		set frame_start [expr $n * $frame]
		set frame_end [expr $frame_start + $n - 1]
		set frame_data [lrange $all_frames $frame_start $frame_end]
		$dsel setbonds $frame_data
		animate goto $frame
		set frame_name [format %05i $frame]
		render snapshot "frame$frame_name.png"
	}
}
