from __future__ import annotations

import math
from typing import TextIO

import numpy as np
import openmm as mm

from ..__rust import Fragment
from ..__simulation import Reporter, Simulation


class LocalMinimizer(Reporter):
    def on_simulation_start(self, simulation: Simulation, continue_sim: bool) -> None:
        pass

    def on_simulation_finish(self, simulation: Simulation) -> None:
        pass

    def __init__(
        self,
        minimizer: LocalGradientDescent,
        minimization_steps: int = 500,
        r_movable: float = 1.0,
        whole_molecule: bool = True,
        harmonic_constraints: bool = True,
        report_every: int = 0,
    ) -> None:
        """
        Local Minimizer. Uses the Reporter API to locally minimize the energy after the modification algorithm runs.

        :param minimizer: The minimization algorithm choice. LocalGradientDescent is currently the only one.
        :param minimization_steps: The number of steps to run the minimization algorithm for.
        :param r_movable: Distance cutoff of how far from reacting atoms counts as local, in nanometers. Set to 0. to only consider the bond graph.
        :param whole_molecule: If True, the entire molecule will be considered local and included in the minimization. If false, only reactant atoms and their bond neighbors will be included.
        :param harmonic_constraints: If True, the constraints will be converted to stiff harmonic bonds during minimization.
        :param report_every: If set to an integer larger than 0, minimization progress details will be written to log files.
        """
        self.minimizer = minimizer
        assert minimization_steps > 0, "Must specify minimization_steps > 0"
        self.minimization_steps = minimization_steps
        self.r_movable = r_movable
        self.whole_molecule = whole_molecule
        self.harmonic_constraints = harmonic_constraints
        self.report_every = report_every
        self.integrator_index = None

    def pre_simulation_start(self, simulation: Simulation) -> None:
        assert simulation.integrator is not None
        self.integrator_index = simulation.integrator.addIntegrator(
            self.minimizer.integrator
        )

    def reset(self, shape: tuple[int, int]) -> None:
        self.minimizer._reset(shape)

    def report(self, handle: TextIO, rem: int) -> None:
        """Write optimization progress to the file handle."""
        handle.write("==============================================\n")
        handle.write(f"remaining steps: {rem}\n")
        handle.write(self.minimizer._report())
        handle.write("\n==============================================\n")

    def on_reaction(
        self, simulation: Simulation, reactions: list[tuple[str, list[Fragment]]]
    ) -> None:
        # save velocities
        assert self.integrator_index is not None
        simulation.info("on_reaction Local Minimizer")
        vel = simulation.context.get_velocities()
        pos, box = simulation.context.get_positions()
        old_integrator = simulation.context.get_current_integrator()
        simulation.context.set_current_integrator(self.integrator_index)
        simulation.top.toggle_softcore(reactions, True)
        if self.harmonic_constraints:
            simulation.system.toggle_constraints_as_harmonic_bonds(True)
        simulation.info("prep done")
        # integrator state setup
        self.reset(vel.shape)
        # set movable
        # atoms in reacting atoms only
        atoms = set()
        for _, frags in reactions:
            for frag in frags:
                for atom in frag.atoms:
                    if atom != -1:
                        atoms.add(atom)
        atoms = simulation.system.populate_neighbors(
            atoms, recursive=self.whole_molecule
        )
        if self.r_movable > 0.0:
            atoms |= box.which_atoms_within_distance(pos, atoms, self.r_movable)

        movable = np.zeros(shape=vel.shape)
        for atom in atoms:
            movable[atom, :] = 1.0

        simulation.info("atoms chosen")
        self.minimizer._set_movable(movable)
        simulation.info("atoms set")

        gcd = math.gcd(self.report_every, 10) if self.report_every > 0 else 10
        remaining = self.minimization_steps
        suffix = f"_minimization{simulation.current_step}.log"
        if self.report_every > 0:
            handle = open(simulation.request_path(suffix), "w")  # noqa: SIM115
            self.report(handle, remaining)

        simulation.info("initial reporting done, minimizing now")
        while remaining > 0:
            c_steps = min(gcd, remaining)
            simulation.context.do_steps(c_steps)
            remaining -= c_steps
            if self.report_every > 0 and remaining % self.report_every == 0:
                self.report(handle, remaining)
            if self.minimizer._has_converged():
                simulation.info("minimizer converged")
                break

        if self.report_every > 0:
            handle.close()

        simulation.info("minimization over")
        simulation.context.set_velocities(vel)
        simulation.info("velocities reset")
        if self.harmonic_constraints:
            simulation.system.toggle_constraints_as_harmonic_bonds(False)
        simulation.top.toggle_softcore(reactions, False)
        simulation.info("preparing to set integrator back")
        simulation.context.set_current_integrator(old_integrator)
        simulation.info("back")


class LocalGradientDescent:
    def __init__(
        self,
        initial_step_size_nm: float = 0.1,
        etol: float = 0.0,
        smoothing_factor: float = 0.1,
    ) -> None:
        """
        Construct a (smoothed) gradient descent minimization integrator.

        :param initial_step_size_nm: Only matters at the start. An adaptive step size is used.
        :param etol: energy tolerance, will stop doing anything once the change in energy reaches this for 1 step. If 0, no convergence check is performed.
        :param smoothing_factor: 0 to 1, the smaller, the smoother but slower convergence, but it could be more stable.

        Modified version, originally from OpenMM Tools
        https://github.com/choderalab/openmmtools/blob/main/openmmtools/integrators.py
        originally distributed under the MIT license,
        Copyright (c) 2015-2019 Chodera lab // Memorial Sloan Kettering Cancer Center.
        """
        assert initial_step_size_nm > 0.0, "initial step size must be larger than 0"
        assert 0.0 < smoothing_factor <= 1.0, "smoothing factor must be between 0 and 1"

        self.global_variables = {
            "step_size": initial_step_size_nm,
            "energy_old": 0,
            "energy_new": 0,
            "delta_energy": 0,
            "accept": 0,
            "fnorm2": 0,
            "eta": smoothing_factor,
            "converged": 0,
            "etol": etol,
            "x_sum": 0,
            "x_sum2": 0,
            "v_sum": 0,
        }

        self.per_dof_variables = {"x_old": 0, "v_old": 0, "est_grad": 0, "movable": 0}

        self.integrator = mm.CustomIntegrator(0.0)

        for k, v in self.global_variables.items():
            self.integrator.addGlobalVariable(k, v)

        for k, v in self.per_dof_variables.items():
            self.integrator.addPerDofVariable(k, v)

        if etol > 0:
            self.integrator.beginIfBlock("converged < 1")

        # Update context state (= Virtual Sites and coupling).
        self.integrator.addUpdateContextState()
        # Constrain positions.
        self.integrator.addConstrainPositions()

        # Store old energy and positions.
        self.integrator.addComputeGlobal("energy_old", "energy")
        self.integrator.addComputePerDof("x_old", "x")
        self.integrator.addComputePerDof("v_old", "v")
        self.integrator.addComputeSum("x_sum", "x")
        self.integrator.addComputeSum("x_sum2", "x*x")
        self.integrator.addComputeSum("v_sum", "v")

        # Take step, re-constraint positions.
        self.integrator.addComputePerDof("est_grad", "(1-eta)*est_grad + eta*f")
        self.integrator.addComputeSum("fnorm2", "est_grad^2")
        self.integrator.addComputePerDof(
            "x", "x+movable*step_size*est_grad/sqrt(fnorm2 + delta(fnorm2))"
        )
        # x was changed, we need to do this again to actually be able to
        # see the energy
        self.integrator.addUpdateContextState()
        self.integrator.addConstrainPositions()

        # Ensure we only keep steps that go downhill in energy.
        self.integrator.addComputeGlobal("energy_new", "energy")
        self.integrator.addComputeGlobal("delta_energy", "energy_new-energy_old")

        # Accept also checks for NaN
        self.integrator.addComputeGlobal(
            "accept", "step(-delta_energy) * delta(energy - energy_new)"
        )

        self.integrator.beginIfBlock("accept = 0")
        # revert DoF positions
        self.integrator.addComputePerDof("x", "x_old")
        self.integrator.endBlock()
        self.integrator.addComputePerDof("v", "v_old")
        # self.addComputePerDof("x", "accept*x + (1-accept)*x_old")

        # Update step size.
        self.integrator.addComputeGlobal(
            "step_size", "step_size * (2.0*accept + 0.5*(1-accept))"
        )

        if etol > 0:
            # check convergence - must be not NaN and delta_energy < 0 and delta_energy > -etol
            self.integrator.addComputeGlobal(
                "converged",
                "delta(energy-energy_new) * step(-delta_energy) * step(delta_energy + etol)",
            )
            self.integrator.endBlock()

    def _set_movable(self, movable: np.ndarray) -> None:
        """Set the degrees of freedom that are movable."""
        self.integrator.setPerDofVariableByName("movable", movable)

    def _has_converged(self) -> bool:
        """Has the minimization converged."""
        return self.integrator.getGlobalVariableByName("converged") == 1

    def _reset(self, shape: tuple[int, int]) -> None:
        """Reset the integrator to its initial state."""
        for k, v in self.global_variables.items():
            self.integrator.setGlobalVariableByName(k, v)

        for k, v in self.per_dof_variables.items():
            assert v == 0
            vals = np.zeros(shape)
            self.integrator.setPerDofVariableByName(k, vals)

    def _report(self) -> str:
        """Return a string that would be reported."""
        return (
            "\n".join(
                f"{k}: {self.integrator.getGlobalVariableByName(k)}"
                for k in self.global_variables
            )
            + "\n"
        )
