source ../../daemon.tcl
daemon_open out.gro out.xtc
daemon_bonds out.bonds

mol modstyle 0 0 CPK 1.5 1.2 12.0 12.0
mol modselect 0 0 resname BDT
mol modcolor 0 0 ColorID 17

mol addrep 0
mol modselect 1 0 name BS1 BS2
mol modcolor 1 0 ColorID 24
mol modstyle 1 0 CPK 1.55 1.25 12.0 12.0

mol addrep 0
mol modselect 2 0 name BS1 BS2
mol modcolor 2 0 ColorID 19
mol modstyle 2 0 VDW 0.4 12.0

display depthcue off
display projection Orthographic
display rendermode GLSL
axes location Off
