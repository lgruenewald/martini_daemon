# Advanced visualization script for martini_daemon

# decodes a int32 list that was written to a file, returns a Tcl list of ints
proc decode_list {filename {compressed 1}} {
	set fd [open $filename r]
	fconfigure $fd -translation binary
	set binary [read $fd]
	close $fd
	if {$compressed} {
		set binary [zlib decompress $binary]
	}
	binary scan $binary "i*" result
	return $result
}

# decodes any tcl obj by a string representation (compressed)
proc decode_compressed {filename} {
	set fd [open $filename r]
	fconfigure $fd -translation binary
	set binary [read $fd]
	close $fd
	set result [zlib decompress $binary]
	return $result
}

# opens a .gro and .xtc with CPK that can be used for visualizing e.g. the bonds
# removes the extra frame from the .gro
# supply molid if opening as a different molid, though the supplied molid
# should be manually verified to be the next one to be added
proc daemon_open {gro xtc {molid 0}} {
	mol new $gro
	mol addfile $xtc waitfor all
	animate delete  beg 0 end 0 skip 0 $molid
}

# colors molecule molid according to the daemon reporter output file npy
proc daemon_color {npy {molid 0} {compressed 1}} {
	
	set asel [ atomselect $molid all ]
	set n [ $asel num ]
	set all_frames [ decode_list $npy $compressed ]
	set len [ llength $all_frames ]
	# number of atoms, assumed to be costant
	set n_frames [expr $len / $n]
	
	for {set frame 0} {$frame < $n_frames} {incr frame 1} {
		# for every frame...
		$asel frame $frame
		$asel update
		set frame_start [expr $n * $frame]
		set frame_end [expr $frame_start + $n - 1]
		# lrange is inclusive on both ends
		set frame_data [ lrange $all_frames $frame_start $frame_end]
		# set custom data
		$asel set user $frame_data
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
		mol selupdate $i $molid 1
	}

	# display some diagnostic data
	puts "daemon_color: loaded $len integers; atom n $n; frame n $n_frames;"
	
}

proc adj_bonds { name element op } {
	# args are vmd_frame molid write, not relevant for us
	# get current frame
	global all_bonds

	set asel [atomselect 0 all]
	set n [$asel num]
	set frame [molinfo 0 get frame]

	# calculate where in the all_bonds we wanna be
	set frame_start [expr $n * $frame]
	set frame_end [expr $frame_start + $n - 1]
	set frame_data [lrange $all_bonds $frame_start $frame_end]

	# update $asel to the current frame and set bond data
	$asel frame $frame
	$asel update
	$asel setbonds $frame_data
}

# displays bonds in frame $frame for molecule $molid according to the
# daemon bonds reporter output file $npy
proc daemon_bonds {npy {molid 0}} {
	global vmd_frame
	global all_bonds
	set all_bonds [ decode_compressed $npy ]
	set llen [ llength $all_bonds ]
	set clen [ string length $all_bonds ]
	puts "daemon_bonds all_bonds: $llen lists $clen chars."

	trace add variable vmd_frame write adj_bonds
	puts "daemon_bonds: vmd_frame trace added."

	adj_bonds a a a
}

