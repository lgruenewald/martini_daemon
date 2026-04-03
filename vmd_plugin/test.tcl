mol new out.gro
mol addfile out.xtc waitfor all
animate delete beg 0 end 0 skip 0 0

load ./toptraj.so
load_toptraj out.toptraj 0 1
