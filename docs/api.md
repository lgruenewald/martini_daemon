This document briefly explains the python API of martini_daemon.
This document is only brief for now, because a lot of these details are up
to change at the moment.

# Importing martini_daemon

The recommended way to import martini_daemon right now is:
```
from martini_daemon.simulation import DaemonSimulation
```

This imports the main class to use.

# Constructing the DaemonSimulation object

The constructor takes two mandatory positional arguments:

- the path of the Topology (.top) file, must be a string
- the path of the Coordinate (.gro) file, must be a string

The constructor takes the following named arguments:

- `T`, default 300.0, the temperature coupling temperature of the system in 
  kelvin, other things like generating velocities is also based on this value
- `p`, default 1.0, the pressure coupling pressure of the system in bars,
  can be set to None to disable pressure coupling
- `dt`, default 20*femtosecond, the timestep
- `max_steps`, default 100, how many reaction steps to do (the D/M algorithm is
  ran once per reaction step)
- `steps_per_step`, default 5000, how many md steps per reaction step
- `traj_path`, default "traj.xtc", the trajectory file to write, currently
  one frame is written once per reaction step
- `out_path`, default "final.gro", the file to write the coordinates at the
  end of the simulation
- `silent`, default False, if True daemon will not write any files nor to the
  stdout, only useful for benchmarking
- `platform`, default None, if the name of an openmm platform is given
  (as a string), it uses that one, with None it will fall back to whichever
  platform openmm considers default
- `minimize_energy`, default True, if True, an energy minimization will be
  performed before doing any md steps
- `generate_velocities`, default True, if True, velocities will be generated
  at the temperature `T`. There is no way yet to read velocities from a file.
- `remove_com_motion`, default True, if True the center of mass motion
  remover force is added to the system
- `epsilon_r`, default 15.0, the epsilon_r value is used for Coulomb
  interactions
- `nonbonded_cutoff`, default 1.1*nanometer
- `include_dir`, default None, daemon will also look for *.itp's in this
  directory when #include-ing from .top files
- `defines`, default {}, extra #define defines to be set during parsing
- `top_logpath`, default None, T* log path, if set it will dump longs there
  that contain information about the current fragment list and the D/M
  algorithm
- `reporters`, default [], a list of reporters (available ones are in the
  submodule martini_daemon.reporters), which can write extra information
  from simulations to files, for example the list of bonds per xtc frame.

# Running simulations

Once a simulation object is constructed, it has the following methods:

- `step()`, which performs `steps_per_step` md steps followed by 
  one reaction step
- `simulate()`, which performs `max_steps` reaction steps

If not silent, both should print out the current step index and the number
of reactions done by the D/M algorithm live to the terminal.
