source ../../daemon.tcl
daemon_open out.gro out.xtc
wait 5
daemon_color out_types.npy

mol modstyle 0 0 CPK 8.0 12.0 12.0 12.0
mol modmaterial 0 0 Transparent
mol modcolor 0 0 ColorID 10
mol modstyle 1 0 CPK 8.0 12.0 12.0 12.0
mol modcolor 1 0 ColorID 1
display projection Orthographic

animate delete  beg 0 end 1 skip 0 0

