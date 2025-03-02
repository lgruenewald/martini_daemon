source ../../daemon.tcl
daemon_open out.gro out.xtc
daemon_bonds out_bonds.npy

mol modstyle 0 0 CPK 2.8 1.3 12.0 12.0
mol modselect 0 0 resname BDT
mol modcolor 0 0 ColorID 3
mol addrep 0
mol modstyle 1 0 QuickSurf 1.600000 0.500000 1.000000 1.000000
mol modselect 1 0 resname BENZ
mol modmaterial 1 0 Transparent
mol modcolor 1 0 ColorID 5
