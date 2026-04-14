"""Context inheriting openmm context"""

import numpy as np
import numpy.typing as npt
import openmm as mm
from openmm.unit import md_unit_system

from ..__rust import PeriodicBox
from .system import System


class Context:
    # automatically read it out from the System and all Force's if there should be a reinitialize!
    # when any state is queried or when there is steps forward, check if reinitializing is needed first, reinitialize
    # and then proceed. MAYBE: have a separate reinitialize only for timing purposes

    def __init__(
        self,
        system: System,
        integrator: mm.Integrator,
        default_box: PeriodicBox,
        platform: mm.Platform | None = None,
        params: dict[str, str] | None = None,
    ):
        """Context object. Created by Simulation automatically, based on the provided `.top` and geometry file.

        :param system: Martini Daemon system linked with this context.
        :param integrator: OpenMM integrator linked with this context. Simulation always builds CompoundIntegrators.
        :param default_box: Default periodic box, used at the start of the simulation.
        :param platform: OpenMM Platform linked with this context.
        :param params: OpenMM Parameters linked with this context.
        """
        system._set_default_pbc(default_box)
        if platform is None:
            self.__context = mm.Context(system.get_openmm_system(), integrator)
        elif params is None:
            self.__context = mm.Context(
                system.get_openmm_system(), integrator, platform
            )
        else:
            self.__context = mm.Context(
                system.get_openmm_system(), integrator, platform, params
            )
        self.system = system
        self.__integrator: mm.Integrator = integrator
        system._bind_context(self)
        self.N = system.num_atoms()

        # before doing any steps or exporting a state, if this is True, reinitialize before proceeding
        self.reinitialize = False

    def do_steps(self, n: int) -> None:
        """Perform integration steps.

        Note: reinitializes if needed automatically.

        :param n: number of integration steps to perform.
        """
        self.__reinitialize()
        if n > 0:
            self.__integrator.step(n)

    def __reinitialize(self) -> None:
        """Automatically called before stepping, reading energies or forces.

        If self.reinitialize was "flagged" to True, rebuilds system and then reinitializes the context.
        """
        if not self.reinitialize:
            return
        self.system._rebuild()
        self.__context.reinitialize(preserveState=True)
        self.reinitialize = False

    def get_current_integrator(self) -> int:
        """If supplied integrator was an OpenMM compound integrator, get the current index.

        Note: Simulation always sets up compound integrators, so generally it will be a compound integrator if using
        the Simulation API.

        :return: Index of current integrator, in the compound integrator.
        """
        assert isinstance(self.__integrator, mm.CompoundIntegrator)
        return self.__integrator.getCurrentIntegrator()

    def set_current_integrator(self, idx: int) -> None:
        """If supplied integrator was an OpenMM compound integrator, set the current index.

        :param idx: Index of integrator within compound integrator to set.
        """
        assert isinstance(self.__integrator, mm.CompoundIntegrator)
        self.__integrator.setCurrentIntegrator(idx)

    def set_positions(
        self, positions: npt.NDArray[np.float64], box: PeriodicBox
    ) -> None:
        """Set atom positions and periodic box vectors in context.

        Constraints, virtual sites, and other interactions that are not periodic in OpenMM are made whole, so that
        their constituting atoms are in a single copy of the periodic space.

        :param positions: Numpy array of positions, in nanometers.
        :param box: PeriodicBox to set.
        """
        if positions.shape != (self.N, 3):
            raise ValueError(
                f"Positions must have shape ({self.N}, 3), got {positions.shape}."
            )
        bonds = self.system.collect_bonds_for_whole()
        bonds.make_whole(box, positions)
        self.__context.setPositions(positions)
        self.__context.setPeriodicBoxVectors(box.a, box.b, box.c)

    def set_velocities(self, velocities: npt.NDArray[np.float64]) -> None:
        """Set atom velocities.

        :param velocities: Numpy array of velocities, in nanometers / picosecond.
        """
        if velocities.shape != (self.N, 3):
            raise ValueError(
                f"Velocities must have shape ({self.N}, 3), got {velocities.shape}."
            )
        self.__context.setVelocities(velocities)

    def generate_velocities(self, temp: float) -> None:
        self.__context.setVelocitiesToTemperature(temp)

    def minimize_energy(self, tolerance: float = 10.0, max_steps: int = 0) -> None:
        """Minimize the potential energy of the system. Calls OpenMM.LocalEnergyMinimizer.minimize().

        :param tolerance: root-mean-square deviation of all forces must be below this value to stop minimization,
            in kJ/mol/nm.
        :param max_steps: maximum number of minimization steps. If 0, minimization will be performed until it
            converges.
        """
        mm.LocalEnergyMinimizer.minimize(self.__context, tolerance, max_steps)

    def apply_constraints(self, tol: float = 1e-10) -> None:
        """Recalculates constraint and virtual site positions. Calls OpenMM.Context.applyConstraints(tol).

        :param tol: distance tolerance for the constraint calculation.
        """
        self.__context.applyConstraints(tol=tol)

    def get_positions(self) -> tuple[np.ndarray, PeriodicBox]:
        """Get current atom positions from the context.

        Positions returned are put in a single copy of the periodic box -- they are not whole!

        Raises an error if positions are beyond what can be stored in the XTC format with default precision.

        :return: positions, in nanometers, and the current periodic box.
        """
        state = self.__context.getState(positions=True)
        vecs = state.getPeriodicBoxVectors()
        box = PeriodicBox(
            [vecs[0].x, vecs[0].y, vecs[0].z],
            [vecs[1].x, vecs[1].y, vecs[1].z],
            [vecs[2].x, vecs[2].y, vecs[2].z],
        )
        pos = state.getPositions(asNumpy=True).value_in_unit_system(
            md_unit_system
        )  # nm
        if np.any(np.abs(pos) > 2147483.0):
            # would be too large to store without remaindering, so likely
            # the system blew up
            raise ValueError(
                "Coordinates too large, your system likely blew up. "
                f"Largest coordinate (abs value) is {np.max(np.abs(pos))}."
            )
        box.move_all_within(pos)
        return pos, box

    def get_velocities(self) -> npt.NDArray[np.float64]:
        """Get current atom velocities from the context.

        :return: velocities, in nanometers/ps.
        """
        state = self.__context.getState(velocities=True)
        return state.getVelocities(asNumpy=True).value_in_unit_system(
            mm.unit.md_unit_system
        )  # nm / ps

    def get_energies(self) -> tuple[float, float, float]:
        """Get the current energies in the system.

        :return: kinetic, potential and total energy of the system, in kJ/mol.
        """
        self.__reinitialize()
        state = self.__context.getState(energy=True)
        pe = state.getPotentialEnergy().value_in_unit_system(
            mm.unit.md_unit_system
        )  # kJ/mol
        ke = state.getKineticEnergy().value_in_unit_system(
            mm.unit.md_unit_system
        )  # kJ/mol
        te = pe + ke
        return ke, pe, te

    def get_forces(self) -> npt.NDArray[np.float64]:
        """Get the current forces acting on each atom.

        :return: forces for each atom, in kJ/(mol nm).
        """
        self.__reinitialize()
        return (
            self.__context.getState(forces=True)
            .getForces(asNumpy=True)
            .value_in_unit_system(mm.unit.md_unit_system)
        )
