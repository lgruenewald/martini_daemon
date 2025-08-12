load ../../tcl/bond_loader.so
source ../../tcl/daemon.tcl
daemon_open out.gro out.xtc
daemon_bonds out.z

animate goto 0
display depthcue off
display projection Orthographic
display rendermode GLSL

