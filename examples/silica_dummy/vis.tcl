load ../../vmd_plugin/toptraj.so

mol new system.gro
mol addfile out.xtc waitfor all
animate delete beg 0 end 0 skip 0 0

load_toptraj out.toptraj 0 1

mol modstyle 0 0 VDW 0.8 12.0
mol modselect 0 0 name SI1
mol modcolor 0 0 ColorID 27
mol addrep 0
mol modselect 1 0 name SV1 SV2 SV3 SV4
mol modstyle 1 0 CPK 1.9 1.3 12.0 12.0
mol modcolor 1 0 ColorID 25
mol addrep 0
mol modselect 2 0 name W
mol modmaterial 2 0 Ghost
mol modstyle 2 0 QuickSurf 1.6 0.5 1.0 1.0
animate goto 0
