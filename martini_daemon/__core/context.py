"""
Context inheriting openmm context
"""
import openmm as mm
import numpy as np
from .system import System
from ..__rust import PeriodicBox
from .utils import collect_bonds_for_whole, make_whole_frame, make_cluster_frame


class Context(mm.Context):
    # automatically read it out from the System and all Force's if there should be a reinitialize!
    # when any state is queried or when there is steps forward, check if reinitializing is needed first, reinitialize
    # and then proceed. MAYBE: have a separate reinitialize only for timing purposes


    def __init__(self, system: System, integrator, default_box: PeriodicBox, platform=None, params=None):
        system._set_default_pbc(default_box)
        if platform is None:
            self.__context = mm.Context(system.get_openmm_system(), integrator)
        elif params is None:
            self.__context = mm.Context(system.get_openmm_system(), integrator, platform)
        else:
            self.__context = mm.Context(system.get_openmm_system(), integrator, platform, params)
        self.system = system
        self.integrator = integrator
        system._bind_context(self)
        self.N = system.atom_count()

        # before doing any steps or exporting a state, if this is True, reinitialize before proceeding
        self.reinitialize = False

    def __reinitialize(self):
        """
        Call this before stepping, reading energies or forces.
        """
        if not self.reinitialize:
            return
        self.system._rebuild()
        self.__context.reinitialize(preserveState=True)

    def set_positions(self, positions: np.ndarray, box: PeriodicBox) -> None:
        if positions.shape != (self.N, 3):
            raise ValueError(f"Positions must have shape ({self.N}, 3), got {positions.shape}.")
        positions = positions.copy()
        bonds = collect_bonds_for_whole(self.system)
        clus = make_cluster_frame(self.N, bonds)
        make_whole_frame(self.N, positions, box, clus)
        self.__context.setPositions(positions)

    def set_velocities(self, velocities: np.ndarray) -> None:
        if velocities.shape != (self.N, 3):
            raise ValueError(f"Velocities must have shape ({self.N}, 3), got {velocities.shape}.")
        self.__context.setVelocities(velocities)

    def generate_velocities(self, temp: float) -> None:
        self.__context.setVelocitiesToTemperature(temp)

    def minimize_energy(self, tolerance=10, max_steps=0) -> None:
        mm.LocalEnergyMinimizer.minimize(self.__context, tolerance, max_steps)

    def apply_constraints(self):
        pos_before = self.__context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(mm.unit.nanometer)
        self.__context.applyConstraints(tol=1e-10)
        pos_after = self.__context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(mm.unit.nanometer)
        return pos_before, pos_after


    def get_positions(self) -> tuple[np.ndarray, PeriodicBox]:
        state = self.__context.getState(positions=True)
        vecs = state.getPeriodicBoxVectors()
        box = PeriodicBox(
            [vecs[0].x, vecs[0].y, vecs[0].z],
            [vecs[1].x, vecs[1].y, vecs[1].z],
            [vecs[2].x, vecs[2].y, vecs[2].z],
        )
        pos = state.getPositions(asNumpy=True).value_in_unit(mm.unit.nanometer)
        if np.any(np.abs(pos) > 2147483.0):
            # would be too large to store without remaindering, so likely
            # the system blew up
            raise ValueError (
                "Coordinates too large, your system likely blew up. "
                f"Largest coordinate (abs value) is {np.max(np.abs(pos))}."
            )
        box.move_all_within(pos)
        return pos, box

    def get_velocities(self) -> np.ndarray:
        state = self.__context.getState(velocities=True)
        return state.getVelocities(asNumpy=True).value_in_unit_system(mm.unit.md_unit_system)  # nm / ps

    def get_energies(self) -> tuple[float, float, float]:
        """
        Kinetic, Potential, Total Energy, in MD units
        """
        self.__reinitialize()
        state = self.__context.getState(energy=True)
        pe = state.getPotentialEnergy().value_in_unit_system(mm.unit.md_unit_system)
        ke = state.getKineticEnergy().value_in_unit_system(mm.unit.md_unit_system)
        te = pe + ke
        return ke, pe, te

    def get_forces(self) -> np.ndarray:
        self.__reinitialize()
        return self.__context.getState(forces=True).getForces(asNumpy=True).value_in_unit_system(mm.unit.md_unit_system)

