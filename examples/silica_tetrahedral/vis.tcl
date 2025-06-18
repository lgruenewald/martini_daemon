load ../../tcl/bond_loader.so
source ../../tcl/daemon.tcl
daemon_open out.gro out.xtc
daemon_bonds out.z

color Display Background white
mol modstyle 0 0 CPK 1.5 1.3 12.0 12.0
mol modselect 0 0 not resname W
mol modcolor 0 0 ColorID 27
mol addrep 0
mol modselect 1 0 name W
mol modmaterial 1 0 Ghost
mol modstyle 1 0 QuickSurf 2.0 0.5 1.0 1.0
animate goto 0
display depthcue off
display projection Orthographic
display rendermode GLSL
