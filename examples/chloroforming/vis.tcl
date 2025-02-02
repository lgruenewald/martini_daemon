source ../../daemon.tcl
daemon_open out.gro out.xtc
wait 5
daemon_color out_types.npy

mol modstyle 0 0 QuickSurf 1.300000 0.500000 1.000000 1.000000
mol modmaterial 0 0 Transparent
mol modcolor 0 0 ColorID 10
mol modstyle 1 0 QuickSurf 1.300000 0.500000 1.000000 1.000000
mol modcolor 1 0 ColorID 1

