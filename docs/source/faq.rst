
Frequently Asked Questions
==========================

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

Specifying the GPU and CPU cores to use
---------------------------------------

Pass ``context_parameters={"DeviceIndex": "0"}`` to the constructor of ``Simulation`` to choose GPU 0.

Selecting CPU cores for the simulation can be done with the ``taskset`` command. For example,
``taskset -c 0-63:2 python3 run.py`` will limit run.py to even numbered cores 0 to 62.
Typically when using a GPU platform, such as CUDA in OpenMM, virtually all the work is performed on the GPU though,
meaning CPU consumption should never be very high when running simulations.

Using Martini Daemon as a ``.top`` file parser for OpenMM
---------------------------------------------------------

Since Martini Daemon implements a large superset of the ``.top`` format, and allows the user to extend it
with custom bond types (see :doc:`/extending`), you might want to use Martini Daemon even for simulations containing
no reactions. If you would like to use the OpenMM Simulation object and OpenMM Reporters, it is possible to export
the OpenMM System and Topology objects from Martini Daemon. Here is a quick script template doing that:

::

    #!/usr/bin/env python3

    import martini_daemon as daemon
    import openmm as mm
    import openmm.app as mmapp

    daemon_sim = daemon.Simulation(
        "system.top", "system.gro",
        sim_name="out",
        # We won't be using Martini Daemon to simulate, so this does not matter
        md_steps=0,
        reporters=[]
    )

    top = daemon_sim.get_openmm_topology()
    sys = daemon_sim.system.get_openmm_system()

    # Do not use daemon simulation any more if using sys on your own!
    del daemon_sim

    # Feel free to do with top and sys as you would with any OpenMM Topology and System instance
    sim = mmapp.Simulation(
        top, sys, mm.VerletIntegrator(0.02), mm.Platform.getPlatformByName("CUDA")
    )

    # At this point, continue as usual with typical OpenMM workflows
    pos = mmapp.GromacsGroFile("system.gro").getPositions(asNumpy=True)
    sim.context.setPositions(pos)

Additional note: if a ``.gro`` file is not specified, Martini Daemon will still parse the ``.top`` file,
but it will not initialize a context. The system returned by ``get_openmm_system()`` will then not
have the default PBC vectors set to the PBC described in the ``.gro`` file. You will have to do that
yourself in that case.
