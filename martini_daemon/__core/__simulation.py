"""
More elegant simulation object that gets passed around everywhere
EVERY COMPONENT SHOULD GET ACCESS TO LOGGING - main log file, screen logging
CUSTOM LOGGER
Current step, current time, reporters, list of open trajectories, log files, reporter files etc inside it, T*, S*, integrator
Request trajectory
Request gro file
Request log file
Request text file
Request uncompressed binary file
Request compressed binary file
improve the passing around of values, remove all default values except in Simulation
Tests should use Simulation not TopParser - clear internal/external
Plugin architecture + split into plugins (mietini, custom H bond, periodic gaussian)
Components / simulation
Topology meta-info – constraints, vsites, coupling, extra forces
EVENT BASED ARCHITECTURE

"""

class Simulation:
    pass