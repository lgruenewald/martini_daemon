Frequently Asked Questions
==========================

Which Martini .top files are supported?
---------------------------------------

More than other OpenMM implementations of ``.top`` (Gromacs Topology format) parsing to date,
but there are some things still missing.

Unsupported features include:

* nrexcl != 1 (number relative exclusions, used to generate exclusions for second or third indirect neighbors, not used in Martini)
* LJ fudge in ``[defaults]`` (not present in Martini's force field itp)
* ``[bondtypes]``, ``[angletypes]``, ``[dihedraltypes]``, ``[constrainttypes]`` (uncommon in Martini)
* Tabulated bonds (uncommon)
* Things you would define in an ``.mdp`` need defining in terms of OpenMM objects (e.g. coupling, PME), which might
  need some extra work in some cases.

An easy method to check if your system is supported is to try it, and see if the energy and forces are correct.

Specifying the GPU and CPU cores to use
---------------------------------------

Pass ``context_parameters={"DeviceIndex": "0"}`` to the constructor of :doc:`/autoapi/martini_daemon/Simulation` to choose GPU 0.
If passing context parameters, it is also recommended to specify the platform (e.g. ``platform="CUDA"``) when
constructing ``Simulation``.

Selecting CPU cores for the simulation can be done with the ``taskset`` command. For example,
``taskset -c 0-63:2 python3 run.py`` will limit run.py to even numbered cores 0 to 62.
Note that when using a GPU platform (e.g. CUDA) in OpenMM,
typically almost all the work is performed on the GPU,
meaning CPU consumption is typically not high when running simulations.

Can I use Martini Daemon as only a ``.top`` file parser for OpenMM?
-------------------------------------------------------------------

Since Martini Daemon implements a large superset of the ``.top`` format, and allows the user to extend it
with custom bond types (see :doc:`/extending`),
there is reasons to try using Martini Daemon even for simulations containing no reactions.
If you would like to use the OpenMM Simulation object and OpenMM Reporters, it is possible to export
the OpenMM System and Topology objects. Here is a quick script template doing that:

::

    #!/usr/bin/env python3

    import martini_daemon
    import openmm as mm
    import openmm.app as mmapp

    daemon_sim = martini_daemon.Simulation(
        "system.top", "system.gro",
        sim_name="out",
        # We won't be using Martini Daemon to simulate, so this does not matter
        md_steps=0,
        reporters=[]
    )

    top = daemon_sim.get_openmm_topology()
    sys = daemon_sim.system.get_openmm_system()

    # Feel free to do with top and sys as you would with any OpenMM Topology and System instance
    # Note: do not use daemon_sim any more, once putting sys in an OpenMM simulation.
    sim = mmapp.Simulation(
        top, sys, mm.VerletIntegrator(0.02), mm.Platform.getPlatformByName("CUDA")
    )

    # At this point, continue as usual with typical OpenMM workflows
    pos = mmapp.GromacsGroFile("system.gro").getPositions(asNumpy=True)
    sim.context.setPositions(pos)

Additional note: if a ``.gro`` file is not specified, Martini Daemon will still parse the ``.top`` file,
but it will not initialize a context. The system returned by ``get_openmm_system()`` will then not
have the default periodic box (PBC) vectors set to the
PBC described in the ``.gro`` file. You will have to do that
yourself in that case.

Periodic Boundary Conditions
----------------------------

The OpenMM FAQ describes how OpenMM handles periodic boundary conditions
at https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions.
This is not quite true for Martini Daemon simulations. This is how
Martini Daemon handles the periodic boundary condition:

* Internally, OpenMM lets things “flow” out of the periodic box, but this is not a problem, the thin wrapper on top of
  OpenMM hides this detail.

    * When a new bond is formed during a reaction, it cannot be guaranteed that the fragments that form that bond
      are in the same copy of the periodic box.
      Bonds, angles, etc. all have the usesPeriodicBoundaryCondition set to True, so this is not a problem.

    * :doc:`/autoapi/martini_daemon/Context` has methods get_position, which puts the particles back into a single
      copy of the periodic box, before returning the positions. This is the method that the Trajectory reporter
      also uses, so trajectories will also have positions within a single copy of the periodic box.

    * Constraints and virtual sites are not PBC aware in OpenMM. This is handled by Martini Daemon by making these
      constructs whole first. This is implemented in :doc:`/autoapi/martini_daemon/Context` method set_positions,
      so this is well encapsulated. Therefore, input geometries with constraints or virtual sites broken across the
      periodic box are supported. Constraints and virtual sites, however, are not allowed
      to be modified during reactions.


Current Limitations
-------------------

* The number of particles cannot change during reactions.

    * If a reaction preserves each bead, it is trivial to keep the old position and velocity, and possibly just do
      a minimization. New beads would need positions and velocities to be defined, which adds complexity.

    * Existing output formats, such as ``.xtc``, and analysis libraries
      are ill-suited for a variable amount of particles, adding more burden to re-create those tools.

    * While adding new particles is more likely in a coarse grained setting, it still is not a common requirement,
      due to the conservation of mass.

    * This limitation can be bypassed with some creativity, such as converting to/from solvent molecules.
      Note: Beads with no interaction with anything else in the system do not currently work.

* Constraints and virtual sites cannot be created or removed during reactions.

    * Would need a new layer of bookkeeping that was not implemented yet, due to internal indices of them changing
      each time one is removed.

    * Constraints and virtual sites are usually used in Martini for strong interactions, such as ring structures,
      which are unlikely to be changed in reactions, so this was not a priority yet.

* Center of mass virtual sites will not change their parameters if the mass of its constructing particles changes
  during a reaction. Center of mass virtual sites that are constructed from other virtual sites (with mass 0) are
  also not supported.

* Pairs, CMAPs will not change parameters automatically if the type of constructing particles changes.
  This can be done manually by removing and re-adding it.

* The charge and mass of beads, if implied from atom type, will not automatically change if the type changes.
  This should be done manually, if it is desired.

* Position restraints are implemented, however they might only work well for orthogonal boxes that do not shrink/grow.

Duplicate exclusion error
-------------------------

OpenMM does not allow the addition of multiple exclusions between two particles. Here is how Martini Daemon handles
this, and when you will get a duplicate exclusion error.

* Exclusions are of course generated also for bonds, based on the nrexcl entry (comes after the molecule name, always 1 for Martini) in ``[moleculetype]``.

* When starting a simulation for the first time, if multiple exclusions are specified in the ``[moleculetype]``
  directive, only one will be added between any pair of particles.

* When a reaction template is executed, which instructs Martini Daemon to create new exclusions, only one new
  exclusion will be added to any pair of atoms which have exclusions defined for them in a reaction template.

* If a reaction template instructs two particles to be excluded, that are already excluded, Martini Daemon will
  not stop it, and OpenMM will raise a duplicate exclusion error. Since this represents an operation which can
  be seen as a common bug in reaction templates (adding a new bond between two beads already bonded), this is
  currently left as is, until a better design is created.

NaN exceptions
--------------

Also see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan.

If this happened right after a reaction,
the following things can be done:

* Check if the distance and angle reaction conditions are sufficiently strict.

    * When forming a bond, it may be necessary to only allow a reaction to happen if the pre-reaction geometry is not
      too far from the equilibrium bond length.

    * Similar reaction conditions may be needed for angles.

* If this fails, try :doc:`/autoapi/martini_daemon/LocalMinimizer`, it may help in some cases.


Are atoms 0-indexed or 1-indexed?
---------------------------------

Both. Depends on where:

* GROMACS file formats, such as itp or gro are 1-indexed.

* In .rx files, to mirror itp, reactants are 1-indexed.

* Internally in OpenMM, everything is 0-indexed.

* Atom indices through Martini Daemon's Python API, and in Martini Daemon-specific output files are 0-indexed.

* Error messages during parsing of itp files are 1-indexed (open an issue if you find an exception to this).

* Some error messages (e.g. in test suite) explicitly state indexing.

When are trajectory frames written?
-----------------------------------

Also see :doc:`/user_guide` section on Reporters.

Trajectory frames are synchronized for all Reporters using the ``on_trajectory_frame`` hook.
They are written:

* First frame before any MD steps were done, unless loading from a checkpoint,

* After each traj_frequency MD steps,

* At the end of simulations (no duplicate if it's exactly at traj_frequency).

``Simulation.trajectory_frame``` represents the current frame index being written,
if in a reporter currently writing a frame.
If a frame is not currently being written, it can be thought of as the number of frames completed.
