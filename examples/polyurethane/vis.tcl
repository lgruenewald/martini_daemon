load ../../tcl/bond_loader.so
source ../../tcl/daemon.tcl
daemon_open out.gro out.xtc
daemon_bonds out.bonds

mol modstyle 0 0 CPK 0.8 1.3 12.0 12.0
mol modselect 0 0 resname DES
mol modcolor 0 0 Name
mol addrep 0
mol modstyle 1 0 CPK 0.8 1.3 12.0 12.0
mol modselect 1 0 resname PEG
mol modcolor 1 0 ColorID 5
mol addrep 0
mol modstyle 2 0 CPK 0.8 1.3 12.0 12.0
mol modselect 2 0 name COH or name NCO
mol modcolor 2 0 ColorID 1
mol addrep 0
mol modstyle 3 0 QuickSurf 1.600000 0.500000 1.000000 1.000000
mol modselect 3 0 resname W
mol modmaterial 3 0 Ghost
animate goto 0
display depthcue off
display projection Orthographic
display rendermode GLSL

