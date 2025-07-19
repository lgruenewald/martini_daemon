# API

The main API entry point of Martini Daemon is the submodule `simulation`. It should be imported as

```
from martini_daemon.simulation import Simulation
```

Many of the simulation parameters usually passed in `.mdp` options should be given to its constructor. It's assumed here that openmm is also imported as

```
import openmm as mm
```

The following arguments can be provided:

- `top_path` (string) - path to the topology `.top` file
- `geom_path` (string), alias `gro_path` - path to the geometry `.gro` file providing the starting coordinates (and optionally velocities)
- `md_steps` (int) - the number of MD steps to perform in total when the `simulate` method is called
- `dm_frequency` (int) - how often should the Detection/Modification (D/M) algorithm be ran, in number of MD steps between each run of the algorithm
- `traj_frequency` (int), alias `xtc_frequency` - how often should the `.xtc` trajectory file (and other reporters which use it) be written to, in number of MD steps between each trajectory frame
- `sim_name` (string) - a prefix given to all output files from the simulation. If present, they will be backed up similar to Gromacs.
- `chk_path` (string) - if given, the checkpoint file at this path will be loaded. Mutually exclusive with `top_path` and `geom_path`.
- `restraint_coord_path` (string) - if given, position restraints will be based on these coordinates. If not given, position restraints will be based on `geom_path`.
- `integrator` (mm.Integrator) - an OpenMM integrator to use. By default a LangevinMiddleIntegrator set to 300 K, 1/ps and 0.02 ps. For options see http://docs.openmm.org/8.2.0/api-python/library.html#integrators.
- `coupling` (list of mm.Force) - list of OpenMM forces, meant for temperature and pressure coupling. By default contains a MonteCarloBarostat at 1 bar and 300 K. For options see http://docs.openmm.org/8.2.0/api-python/library.html#forces.
- `remove_com_motion` (bool) - whether center of mass motion should be removed. Default: true. Please do not remove com motion through the coupling option, as this sets a flag in S* and affects how VariablesReporter reports temperature by changing the number of degrees of freedom.
- `platform` (None | string | mm.Platform) - if None, the default OpenMM platform is chosen by the library. An OpenMM platform object can be given, or alternatively the name of the platform as a string. Available platforms are "Reference", "CPU", "OpenCL", "CUDA", "HIP". See http://docs.openmm.org/8.2.0/userguide/library/01_introduction.html#platforms for more info.
- `context_parameters` (dict) - The parameters passed to the OpenMM platform. Common use case would be to limit which GPU OpenMM should run on, which can be done by passing `{"DeviceIndex":"0"}` to e.g. run exclusively on GPU 0. See http://docs.openmm.org/8.2.0/userguide/library/04_platform_specifics.html for specifics. Must specify platform to specify context parameters.
- `nonbonded_force` (Nonbonded) - The Martini Daemon force to use as the nonbonded force. Can be specified to e.g. change the nonbonded force cutoff or epsilon. Import `from martini_daemon.forces.nonbonded import NonBonded` and pass `NonBonded(epsilon_r=custom_epsilon_value, cutoff_nm=custom_cutoff)` to set a custom epsilon and cutoff.
- `include_dir` (string) - if specified, #include in the `.top` file will also search in this directory. If not specified, the environment variables GMXDATA and GMXBIN are searched to find the `gromacs/top` folder in the active Gromacs installation and use that as the default. If that fails, it defaults to `/usr/local/gromacs/share/gromacs/top`.
- `defines` (dict[string, string]) - a list of #define key-value pairs for pre-processing `.top` files. By default the key DAEMON is set to "1". Both keys and values should be strings. Set things to "1" to pass #ifdef. Example: {"FLEXIBLE": "1"} to use the flexible versions of constraints found in many parameter files.
- `reporters` (list[Reporter]) - a list of reporters which are called at points of simulations to write information to output files. Example reporters can be found in `martini_daemon.reporters`. See the reporter section of this document for more examples. Note, OpenMM reporters don't work out of the box, Reporter is a Martini Daemon type.
- `neighbor_cutoff` (float) - defaults to 1.1. The cutoff for the detection algorithm neighbor search. Should be increased if r_max for a multi molecular reaction is approaching 1.1.
- `force_reinitialize` (bool) - only used for benchmarking, don't use in production simulations. If set to true, the OpenMM context is reinitialized every step.
- `max_absolute_rate` (float) - if set, relative rate controlled reactions relative rate will be limited to this value. See the section on rate control for more details.
- `rate_highest_probability` (float) - between 0.0 and 1.0, default 1.0. All rate controlled reactions will have the highest probability of this to succeed. See the section on rate control for more details.

Once the Simulation class is constructed, the following methods can be called:

- `minimize_energy()`
- `generate_velocities(T)`
- `step(n)` - performs n steps
- `simulate()` - performs as many steps as specified during the construction of the object, running the Detection/Modification and writing output as appropriate.

Martini Daemon uses OpenMM in the background, however some of the
internal details are hidden. This is because of the bookkeeping
required to support reactions. The bookkeeping is constructed
during the parsing of the system. Martini Daemon also keeps the
state of the OpenMM simulation up to date with the bookkeeping
at all times. This currently means, that only the Martini Daemon
API should be used. If adding further forces, Martini Daemon will
not be aware of them.

# Graphs

- graph input description

- rough description of the algorithm
  - eagerness
  - only useful information for users

## Interaction filters

- `bond` - any bond listed here
- `harmonic_bond` - bond type 1 and 6
- `g96_bond` - bond type 2
- `morse_bond` - bond type 3
- `cubic_bond` - bond type 4
- `connection` - bond type 5
- `fene_bond` - bond type 7
- `distance_restraint` - bond type 10
- `constraint` - constraint type 1 and 2

- `angle` - any angle listed here
- `harmonic_angle` - angle type 1
- `g96_angle` - angle type 2
- `cross_bond_bond` - angle type 3
- `cross_bond_angle` - angle type 4
- `urey_bradley` - angle type 5
- `quartic_angle` - angle type 6
- `linear_angle` - angle type 9
- `restricted_angle` - angle type 10

- `dihedral` - any dihedral listed here
- `proper_dihedral` - dihedral type 1, 4 and 9
- `improper_dihedral` - dihedral type 2
- `rb_torsion` - dihedral type 3 and 5
- `restricted_dihedral` - dihedral type 10
- `combined_bending_torsion` - dihedral type 11

- `virtual_site` - any virtual site, matches the virtual particle as well as constructing particles
- `vsite` - alias for virtual_site
- `vsite1` - virtual_sites1 type 1
- `vsite2` - virtual_sites2 type 1
- `2fd` - virtual_sites2 type 2
- `vsite3` - virtual_sites3 type 1
- `3fd` - virtual_sites3 type 2
- `3fad` - virtual_sites3 type 3
- `3out` - virtual_sites3 type 4
- `4fdn` - virtual_sites4 type 2
- `com` - virtual_sitesn type 2
- `weighed_average` - virtual_sitesn type 1 and 3

- `pair` - only matches pairs
- `cmap` - only matched cmap

# Reaction templates

- describe all options and syntax

# Rate control prototype

## Simulation parameters

- per simulation:
  - `max_absolute_rate` - optional, `absolute_rate` gets capped at this value, in units of reactions per D/M step
  - `highest_probability` - mandatory with default value 1 currently, a number between 0 and 1, the probability of a reaction happening for every reaction will get multiplied by this number, so it can be used to scale all reactions.
- per reaction:
  - `relative_rate` - rate constant `k` for the reaction, relative to `absolute_rate`, therefore contains no units of time, for bimolecular or larger order reactions: concentration units of molecules per simulation box

## State kept

- per simulation:
  - `absolute_rate` - multiply this by each `relative_rate` to get the `k` rate constant that the algorithm aims for
- per reaction:
  - `observed_rate` - connected to the frequency of passing geometry conditions of a reaction

## Goals and limitations

- Non-goal: realistic absolute rates
- Goal: correct instantaneous relative rates between the reactions in the system at a specific point of the simulation
- Built on the assumption, that the system contains all reactants evenly mixed, and that all matches of a single class of reactant (a single graph) react at roughly similar rates
- Limitations: `absolute_rate` can differ during the reaction, so the rate of a reaction can change (even by orders of magnitude) during the same simulation, for example if new, much much slower reactions start happening during the simulation. Generally this rate limiting is most useful if the nature of the system (what types of reactions happen) doesn't change that much during a simulation from the initial point. Or in some cases, it can mean that the initial rates will be different from most of the simulation (e.g. dimerization and breakage - before there is anything to break, we can't correctly set the absolute rate of dimerization to get the good relative rate with breakage, if breakage is much slower).
- Limitation: relative rates are meaningless with only a single reaction. Beware of unimolecular reactions with only a rate condition - if there is only a single breaking reaction in the system it will proceed at 100% speed (all done in a single D/M run) if `max_absolute_rate` or `highest_probability` are not set to limit it.
- Note: Uni and bimolecular reaction actual relative rates to eachother will be concentration dependent, since a dimensional analysis will reveal that `absolute_rate` has a unit problem
  - but this is expected, since saying that reaction A (unimolecular) should be 50x faster than reaction B (bimolecular) at all reactant concentrations is impossible
- currently all reactions in a simulation get tied to the same `absolute_rate`, it's not possible to "group them" (e.g. reactions A and B have a specific relative rate to eachother, independent from reactions C and D), this could be lifted later, maybe?

## Algorithm

- T* contains a value `absolute_rate`, which gives the dimension
  of rate constants (to get a real `k` for every reaction)
- initially (first D/M run) it is not set, the first D/M run is considered a warm up and no rate controlled reactions will happen
- After each D/M run, `observed_rate` is calculated:
  - as reactions during the last step, divided by `relative_rate` and the product of all concentrations (in units of molecules per simulation box) - this way it is converted into the units of `absolute_rate` and already scaled by the `relative_rate`!
  - the `observed_rate` value is kept smoothed using an exponential smoothing (very simple smoothing formula that smoothes the values and slopes for values)
  - if there are no reactants for the reaction, its `observed_rate` is not updated
  - initially `observed_rate` is not set
- `absolute_rate` is set to the smallest of `observed_rate` * `highest_probability`, or `max_absolute_rate` (if specified, and if it is smaller than all of them)
- in the detection algorithm, if a reaction is rate controlled but `absolute_rate` or the reaction specific `observed_rate` are not set, the reaction is rejected
- in the detection algorithm, the probability of accepting the reaction is `absolute_rate` / `observed_rate`
  - Note: for the slowest reaction, `absolute_rate` = `observed_rate`, therefore this probability will be 100% (unless `max_absolute_rate` is smaller)
  - for every other reaction, `observed_rate` is bigger than `absolute_rate`, by a factor that estimates the relative frequency of geometry conditions becoming true between the two reactions

# Reporters

- for each:
  - API + if helper present, its API
  - output format

# Helpers

- for each:
  - API
  - explanation

# OpenMM and Martini Daemon

- currently implemented subset of Martini

## Periodic boundary condition

The OpenMM FAQ describes how OpenMM handles periodic boundary conditions at https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions. This is not quite true for Martini Daemon simulations.

TODO

# Limitations

The number of particles cannot be changed during reactions.

Constraints and virtual sites cannot be created or removed during
reactions. Cannot be created because of periodic boundary conditions
inside openmm. In theory to support adding, the molecules would need
to be made whole across the pbc first. To support removal, the
indexing would need to be kept track of too. Currently the parameters
for them also cannot be changed, as parameter changes are done using
removal and readding currently.

Fragments that overlap can't react.

There is no checking for duplicate exclusions created between two
particles during reactions. If a reaction adds an exclusion (or a
bond that excludes) during a reaction between two atoms that are
already excluded (possibly through a bond), OpenMM will raise an
exception. Don't just add bonds to pairs of atoms already bonded.
Exclusions are deduplicated within the same moleculetype or same
reaction, though, so adding multiple bonds at once or excluding a
bonded pair of atoms is fine.

Reaction constraints and product interactions specified by the user
should lead to forces during bond formation that maintain the
numerical stability of the system, and it's the user's responsibility
to ensure this.

Center of Mass virtual sites will not change parameters if the
constructing particle mass changes during a reaction.

Only periodic boxes with 90 degree angles are supported.

# Known bugs

Pairs and CMAP are broken.

# Programmer's guide

## Code overview

- a few sentences about each file or folder

## Adding new forces

## Adding new reaction conditions
