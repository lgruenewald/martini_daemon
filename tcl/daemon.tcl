# Advanced visualization script for martini_daemon

# opens a .gro and .xtc with CPK that can be used for visualizing e.g. the bonds
# removes the extra frame from the .gro
# supply molid if opening as a different molid, though the supplied molid
# should be manually verified to be the next one to be added
proc daemon_open {gro xtc {molid 0}} {
	mol new $gro
	mol addfile $xtc waitfor all
	animate delete  beg 0 end 0 skip 0 $molid
}

# decodes any tcl obj by a string representation (compressed), depends on zlib
proc decode_compressed {filename} {
	set fd [open $filename r]
	fconfigure $fd -translation binary
	set binary [read $fd]
	close $fd
	return [zlib decompress $binary]
}

# tracing function, depends on bonds_frame being loaded from libdaemon.so
proc adj_bonds { name element op } {
	global trace_molid
	set asel [atomselect $trace_molid all]
	set n [$asel num]
	set frame [molinfo $trace_molid get frame]

	# update $asel to the current frame and set bond data
	$asel frame $frame
	$asel update
	$asel setbonds [bonds_frame $frame]
}

# pure tcl adj_bonds
proc adj_bonds_tcl { name element op } {
	# args are vmd_frame molid write, not relevant for us
	# get current frame
	global all_bonds
	global trace_molid

	set asel [atomselect $trace_molid all]
	set n [$asel num]
	set frame [molinfo $trace_molid get frame]

	# calculate where in the all_bonds we wanna be
	set frame_start [expr $n * $frame]
	set frame_end [expr $frame_start + $n - 1]
	set frame_data [lrange $all_bonds $frame_start $frame_end]

	# update $asel to the current frame and set bond data
	$asel frame $frame
	$asel update
	$asel setbonds $frame_data
}

# setup trace to display bonds
proc daemon_bonds { npy {molid 0} } {
	global vmd_frame
	global trace_molid
	global all_bonds

	set libdaemon_loaded [
		expr [llength [info commands bonds_load]] \
		+ [llength [info commands bonds_frame]] >= 2
	]
	set zlib_loaded [
		expr [llength [info commands zlib]] >= 1
	]

	set asel [atomselect $molid all]
	set n_atoms [$asel num]
	set n_frames [molinfo $molid get numframes]
	set trace_molid $molid
	if {$libdaemon_loaded} {
		puts "daemon_bonds: using libdaemon.so version."
		bonds_load $npy $n_atoms $n_frames
		trace add variable vmd_frame write adj_bonds
		puts "daemon_bonds: vmd_frame trace (adj_bonds) added."
		adj_bonds a a a
	} elseif {$zlib_loaded} {
		puts "daemon_bonds: using pure tcl version (with zlib)."
		set all_bonds [ decode_compressed $npy ]
		trace add variable vmd_frame write adj_bonds_tcl
		puts "daemon_bonds: vmd_frame trace (adj_bonds_tcl) added."
		adj_bonds_tcl a a a
	} else {
		error "Error: neither bonds_load nor zlib exist as tcl commands. Please compile and load libdaemon.so or use a vmd installation with zlib."
	}
}
