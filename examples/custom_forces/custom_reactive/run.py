#!/usr/bin/env python3

from martini_daemon import simulation
from martini_daemon.reporters.bond_reporter import BondReporter
import openmm as mm

sim = simulation.Simulation(
    top_path="system.top", gro_path="system.gro",
    sim_name="out",
    reporters=[
        BondReporter(),
    ],
    integrator=mm.LangevinMiddleIntegrator(298 * mm.unit.kelvin, 1 / mm.unit.picosecond, 0.01 * mm.unit.picosecond),
    coupling=[
        mm.MonteCarloBarostat(1*mm.unit.bar, 298 * mm.unit.kelvin)
    ],
    md_steps=1500000, dm_frequency=0,
    xtc_frequency=5000,
    platform="HIP",
    defines={
        "REACT": "1"
    }
)
sim.minimize_energy()
sim.generate_velocities(298)
sim.simulate()
