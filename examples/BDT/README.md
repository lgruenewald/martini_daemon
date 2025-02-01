A system of 100 benzene dithiol molecules in benzene.

Note: the angle cutoff reaction condition is essential for system stability.

Files:

- BDT.itp - itp of BDT monomer and dimer, as parametrized
- BDT.gro, dimer.gro - single monomer and dimer coordinate files
- benzbox.gro - equilbirated box of martini benzene
- **system.gro** - input file for the system
- **system.top** - topology for system.gro
- **dimerize.rx** - dithiol formation and breakage
- **run.py** - script to run daemon on the system

How to run:

1. Install martini_daemon
2. `./run.py` while in this directory

How to view bonds:

1. open VMD
2. source daemon.tcl (root of this repository)
3. load the trajectory and delete the extra frame for the loaded .gro file (the helper proc `daemon_open coords.gro traj.xtc` does this)
4. load the bond information for a specific frame (call the proc `daemon_bonds out_bonds.npy 0` to load the bonds from frame 0 of out_bonds.npy)

