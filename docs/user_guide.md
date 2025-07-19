# Installation

First, it's recommended to install OpenMM (http://docs.openmm.org/latest/userguide/application/01_getting_started.html) with CUDA (if nvidia) or HIP (if AMD) support. A conda environment or a virtual environment is recommended. Second, clone the repository and install using pip:

```
git clone https://github.com/lgruenewald/martini_daemon
cd martini_daemon
pip install .
```

# Running a simulation

Simulations in Martini Daemon are ran using python run scripts.
These contain calls to Martini Daemon's API, specifying simulation
parameters and input files. A simple energy minimization run script
is provided below, which should be adjustable to meet various
needs.

```
from martini_daemon.simulation import Simulation
from martini_daemon.reporters.variables_reporter import VariablesReporter
import openmm as mm

# equilibration
eq = Simulation(
  # input files
  top_path="system.top", geom_path="system.gro",
  # output files will be named after this
  sim_name="eq",
  # total MD steps to take
  md_steps=400000,
  # how often to write to the .xtc
  traj_frequency=5000,
  # no reactions during equilibration
  dm_frequency=0,
  # report thermodynamic variables, other reporters go here too...
  reporters=[
    VariablesReporter()
  ],
  # integrator+temperature coupling in one
  integrator=mm.LangevinMiddleIntegrator(
    # temperature
    298 * mm.unit.kelvin,
    # collision frequency
    1.0 / mm.unit.picosecond,
    # timestep
    0.02 * mm.unit.picosecond
  ),
  # coupling
  coupling=[
    mm.MonteCarloBarostat(
      # desired pressure
      1.0 * mm.unit.bar,
      # reference temperature
      298 * mm.unit.kelvin
    )
  ],
  # use CUDA with nvidia GPUs
  platform="CUDA"
  # Run on GPU 1 only
  context_parameters={"DeviceIndex": "1"}
)
# minimize energy first, saving the minimized coordinates to min.gro
eq.minimize_energy(out="min.gro")
# generate velocities at 298 K
eq.generate_velocities(298)
# run equilibration
eq.simulate()
```

# Including reactions

1. determine mechanism

2. parametrize reactant and product in mechanism

3. write down difference in terms of new bonds, angles... and broken bonds, angles...

4. make a graph for the reactant

5. check if the graph works

6. make a reaction template

7. run simulation

8. check if reactions work

# Reporting

- Variables Reporter

- Checkpoint Reporter

- Atom Reporter

- Bond Reporter

# Analysis

- Monomer helper

- PBC whole

# Visualization

- VMD helper
