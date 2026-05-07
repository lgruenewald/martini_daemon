Extending Martini Daemon
========================

Martini Daemon is intended to be extensible. Extending is possible in multiple ways:

* Creating new reporters,
* Creating new bond, angle, etc. types,
* Changing the NonBonded force,
* Adding new directives to the ``.top`` parser.

Creating new reporters
----------------------

New reporters should inherit the ``Reporter`` class and override various of its methods. Each of the overridden
methods implements a callback that is called by ``Simulation`` at certain points during Martini Daemon simulations.
:doc:`/autoapi/martini_daemon/Reporter` can be checked for an up-to-date list of possible callbacks. In this section,
a specific example will be shown, along with explanations. This example will be CustomReporter, which will call
methods of an imaginary CustomWriter class. This imaginary class is used to hide the complexity of file format
writers, and to focus on interfacing with Martini Daemon.

::

    from martini_daemon import Reporter, Simulation

    class CustomReporter(Reporter):
        def __init__(self):
            self.writer: CustomWriter | None = None

One of the main methods to implement is ``on_simulation_start``. Reporters during this callback should open
their output files and write their headers. If the ``continue_sim`` argument is true, reporters should open
their output files for appending, and they should truncate the pre-existing output to be aligned with the current
frame, as can be read from Simulation (``Simulation.current_step`` for the current MD step, ``Simulation.time_ps`` for
the current time in ps, ``Simulation.trajectory_frame`` for the number of trajectory frames in the past).
The output file paths should be obtained using ``Simulation.request_path``, since all Reporters should
have the simulation name as a prefix to all output file paths. ``Simulation.request_path`` will also back the file
found at the path up and inform the user of this. If the copy argument of ``Simulation.request_path`` is True,
the backup process is a copy rather than a move. This should be True if ``continue_sim`` is True.

::

        def on_simulation_start(self, simulation: Simulation, continue_sim: bool):
            title = simulation.system.additional_data.get("title")
            if continue_sim:
                self.writer = CustomWriter(
                    path=simulation.request_path(".custom", True),
                    mode="a"
                )
                # Imaginary method that truncates an existing file to a number of frames
                self.writer.truncate(simulation.trajectory_frame)
            else:
                self.writer = CustomWriter(
                    path=simulation.request_path(".custom"),
                    mode="w"
                )
                # Imaginary method that writes a once-per-file header
                self.writer.write_header(title)

The main callback called every trajectory frame is ``on_trajectory_frame``. Generally, this is called at the
start of each simulation, as well as every time after ``traj_frequency`` steps have passed.
Reporters should read out the required data from the ``Simulation`` object, and its fields such as ``Simulation.context``
or ``Simulation.system``. Reporters generally should strive to flush their output after trajectory frames,
so analysis of still running simulations is possible. For our example, the periodic box, positions and list of bonds
is reported at each frame, using CustomWriter's methods.

::

        def on_trajectory_frame(self, simulation):
            pos, box = simulation.context.get_positions()
            # box is a PeriodicBox instance, CustomWriter might expect e.g. the unit cell vectors as a 3x3 numpy array
            box_np = np.array([box.a, box.b, box.c])
            # collect bonds that pass the filter "bond" or "vsite", convert it to a python list of bonds
            bond_list: list[tuple[int, int]] = simulation.system.collect_bonds(["vsite", "bond"]).to_list(),

            assert self.writer is not None
            # Imaginary method that writes frame metadata, atom positions, periodic box and list of bonds
            self.writer.new_frame(
                simulation.trajectory_frame, # the current index of this frame / how many frames were completed before
                simulation.current_step, # current MD step
                simulation.time_ps, # simulation time
                simulation.system.num_atoms(), # number of atoms
                box_np,
                pos,
                bond_list
            )
            # Recommendation: flush every frame
            self.writer.flush()

Resource cleanup happens during ``on_simulation_finish``. It is the assumption that many Reporters will own
open file handles. This gets called when ``Simulation.finish()`` is called either manually, or when it goes
out of scope in a ``with Simulation(...) as sim:`` block. Open file handles owned by the Reporter should be closed
during this method. The simulation object is still available to query, should some information from it be necessary.

::

        def on_simulation_finish(self, simulation):
            assert self.writer is not None
            self.writer.close()

Creating new interactions
-------------------------

For this section, the implementation of Harmonic Bond in Martini Daemon is explained. It can be modified to
get custom bond, angle or dihedral types. Note: this part of the API is newer and still up to change.

::

    import openmm as mm

    from martini_daemon import BondedForce, register_available_force, register_bond_type

    @register_bond_type(type_=1, args=["float", "float"], is_excl=True)
    @register_bond_type(type_=6, args=["float", "float"], is_excl=False)
    @register_available_force
    class HarmonicBond(BondedForce):
        @classmethod
        def _add_to_force(
            cls, force: mm.Force, members: list[int], params: list[float]
        ) -> None:
            assert isinstance(force, mm.HarmonicBondForce)
            force.addBond(*members, *params)

        def _parse(self, members: list[int], params: list[float]) -> list[float]:
            return params

        @classmethod
        def uses_pbc(cls) -> bool:
            return True

        def delta_degrees_of_freedom(self) -> int:
            return 0

        def _set_force_obj(self) -> mm.Force:
            return mm.HarmonicBondForce()

        filters = {"bond"}

        @classmethod
        def get_name(cls) -> str:
            return "harmonic_bond"

Custom bonded interactions should inherit :doc:`/autoapi/martini_daemon/BondedForce`. There is a number of methods
that must be overridden:

* ``register_bond_type`` is a decorator defined for the ``[bonds]`` directive. It binds a specific bond type number to
  a BondedForce subclass, and validates the number and type of arguments specified in ``.itp`` files. Other decorators
  exist for ``[angles]`` and ``[dihedrals]`` directives. As seen in the example, in principle, multiple bond type
  numbers can be associated with a single force. For bonds, ``is_excl`` is also a necessary argument, which specifies
  whether exclusions should be generated by this bond, according to the ``nrexcl`` rule for the moleculetype.

* ``register_available_force`` is a decorator that allows the automatic addition of this force to any Martini Daemon
  System, if a molecule type references it.

* ``_add_to_force`` should call the OpenMM function that adds a new interaction, given members and parsed parameters.
  It can be assumed that the force argument is the same as returned by ``_set_force_obj``, but may be modified.
  This method should have no side effects, other than mutating the passed force,
  as it will be called lazily and possibly multiple times by Martini Daemon.

* ``_parse`` should be used to convert parameters from how they are in input files to how they are stored within this BondedForce instance,
  Some parsing, such as conversion of degrees to radians, or the conversion of ``.itp`` atom indices to OpenMM System indices is already
  performed. Note, that members are only provided as a reference, as in rare cases, querying information from ``self.system`` for
  those members may be necessary.

* ``uses_pbc`` should return whether this force is constructed in a periodic box aware manner (usually True).
  If True, by default, the OpenMM method setUsesPeriodicBoundaryConditions will be called on the force object.

* ``delta_degrees_of_freedom`` should return if the addition of this force changes the number of degrees of freedom in the system (usually 0).

* ``_set_force_obj`` should return the OpenMM force object. This is also called lazily. It may read things out
  from self.system.

* ``filters`` should be a set of filters other than the name, which can be used to reference it. Generally, filters
  such as ``bond``, ``angle`` or ``dihedral`` should be specified, as well as other aliases that a force may be
  known by. These are not unique per force.

* ``get_name`` should return a unique name for the interaction.

Changing the NonBonded force
----------------------------

In broad terms :doc:`/autoapi/martini_daemon/NonBonded`  (found in the repo at ``/martini_daemon/__forces/nonbonded.py``)
works as such:

* The base NonBonded force is a shifted LJ and reaction field electrostatic, with a cutoff. Additionally, it has soft core
  parameters, that can be "toggled on" by changing the values away from the default.
* The parsed nonbonded params directive result from parsing the ``.top`` is put into ``system.additional_data["nb_types"]``.
  This is read by ``NonBonded`` into a Discrete 2D table of OpenMM.
* ``ExclusionHelper`` handles storing exclusions, as well as adding reaction field corrections to the energy.

There is no proper API to only change the NonBonded function. Nevertheless, it can be changed.

* Make a subclass of NonBonded that overrides everything. Copy the implementation and edit it.
* Pass the modified class as the argument ``nonbonded`` to Simulation.

Adding new directives
---------------------

In broad terms:

* Custom directives should inherit the :doc:`/autoapi/martini_daemon/Directive` class.
* Custom directives need to be registered using ``@register_directive``.
* Each directive can specify which other directive should its parent be. The root :doc:`/autoapi/martini_daemon/GromacsTopFile`
  is a likely candidate. The parent directive instance is available under the ``parent`` argument passed to the constructor.
* Each occurrence of a directive during parsing will lead to the instantiation of the registered Directive subclass.
* Temporary state can exist within ``self``.
* The output of parsing is modelled by mutating ``parent``, or its fields (e.g. ``system``).
* A good place to put per-simulation data is ``system.additional_data``, which is a dictionary that stores things
  the core of Martini Daemon need not be aware of, but all the different moving components can still access and mutate
  it. Its lifetime is the current simulation.
