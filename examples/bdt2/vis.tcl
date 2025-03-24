source ../../daemon.tcl
daemon_open out.gro out.xtc
daemon_bonds out_bonds.npy

mol modstyle 0 0 CPK 2.8 1.3 12.0 12.0
mol modselect 0 0 resname BDT
mol modcolor 0 0 ColorID 3
