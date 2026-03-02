import openmm as mm
import numpy as np
from .daemon_integrator import DaemonIntegrator
import math
from ..utils import pdist
from ..reporters.bond_reporter import collect_bonds
from ..helpers.cluster import make_cluster_frame


class LocalMinimizingIntegrator(DaemonIntegrator):
    """
        Inspired by
        https://github.com/choderalab/openmmtools/blob/main/openmmtools/integrators.py
        GradientDescentMinimizationIntegrator
    """

    def __init__(
        self, dt_ps, T_K, friction_ps1, minimizer, minimization_steps=500,
        report_every=0, check_convergence_every=10,
        langevin=True, equilibration_length=0, subdivision=4,
        r_movable = 0., whole_molecule=False, toggle_couplings=False,
        max_retry = 0, harmonic_constraints=False
    ):
        """
        ToggleCouplings is experimental, don't use that one
        """
        self.dt = dt_ps
        self.T = T_K
        self.friction = friction_ps1
        self.integrator = mm.CompoundIntegrator()
        if langevin:
            self.integrator.addIntegrator(
                mm.LangevinMiddleIntegrator(
                    T_K, friction_ps1, dt_ps
                )
            )
        else:
            self.integrator.addIntegrator(
                mm.VerletIntegrator(dt_ps)
            )
        self.minimizer = minimizer
        self.integrator.addIntegrator(minimizer)
        assert minimization_steps > 0, "Must specify minimization_steps > 0"

        if langevin:
            self.integrator.addIntegrator(
                mm.LangevinMiddleIntegrator(
                    T_K, friction_ps1, dt_ps / subdivision
                )
            )
        else:
            self.integrator.addIntegrator(
                mm.VerletIntegrator(dt_ps / subdivision)
            )

        self.minsteps = minimization_steps
        self.report_every = report_every
        self.check_convergence_every = check_convergence_every
        self.remaining_eq_steps = 0
        self.eq_steps = equilibration_length
        self.subdivision = subdivision
        self.r_movable = r_movable
        self.whole_molecule = whole_molecule
        self.toggle_couplings = toggle_couplings
        self.max_retry = max_retry
        self.harmonic_constraints = harmonic_constraints

    def reset(self, shape):
        self.minimizer.reset(shape)

    def report(self, i, rem, system, op="a"):
        with open(f"{self.sim_name}_minimization{i}.log", op) as f:
            f.write("==============================================\n")
            f.write(f"Simulation frame {i}, remaining steps {rem}\n")
            f.write("----------------------------------------------\n")
            f.write(self.minimizer.report())
            f.write("----------------------------------------------\n")
            f.write(system.get_energies_and_forces_by_group())
            f.write("\n==============================================\n")


    def set_reactions(self, reactions, system, top, i, retries_so_far=0):
        # save vels
        state = system._context.getState(
            velocities=True,
            positions=True,
            energy=True
        )
        pos = state.getPositions(
            asNumpy=True
        ).value_in_unit(mm.unit.nanometer)
        vels = state.getVelocities(
            asNumpy=True
        )
        before_energy = state.getPotentialEnergy()
        box = np.array(system.get_box(state))
        self.integrator.setCurrentIntegrator(1)
        if self.harmonic_constraints:
            system.constraint.constraints_to_harmonic_bonds(True)
            system.reinitialize()
        if self.toggle_couplings:
            system.coupling(False)
            system.reinitialize()
        # integrator state setup
        self.reset(vels.shape)
        # set movable
        movable = np.zeros(shape=vels.shape)
        atoms = set()
        for (frags, _) in reactions:
            for frag in frags:
                for atom in frag.atoms:
                    if atom != -1:
                        atoms.add(atom)
        if not self.whole_molecule:
            # old code of 1 depth search
            for atom in atoms.copy():
                for inter in top.interaction_list[atom]:
                    for member in inter.get_members():
                        atoms.add(member)
        else:
            # make entire molecules whole, based on bonds, constraints, virtual sites graph
            _, bonds = collect_bonds(system.len_atoms(), system)
            clus = make_cluster_frame(bonds, len(vels))
            done = set()
            for atom in atoms.copy():
                c = clus[atom]
                assert c > 0
                if c in done:
                    continue
                for natom in np.argwhere(clus == c):
                    assert natom.shape == (1,)
                    atoms.add(natom[0])
                done.add(c)

        for atom in atoms:
            movable[atom, :] = 1.

        if self.r_movable > 0:
            for i, cpos in enumerate(pos):
                for j in atoms:
                    if movable[i, 0] == 1.:
                        break
                    base = pos[j, :]
                    if pdist(cpos, base, box) < self.r_movable:
                        movable[i, :] = 1.

        self.minimizer.setPerDofVariableByName(
            "movable",
            movable
        )

        if self.report_every > 0:
            gcd = math.gcd(self.report_every, self.check_convergence_every)
        else:
            gcd = self.check_convergence_every
        if gcd == 0:
            gcd = self.minsteps

        remaining = self.minsteps

        if self.report_every > 0:
            self.report(i, remaining, system, op="w")

        while remaining > 0:
            csteps = min(gcd, remaining)
            self.integrator.step(csteps)
            remaining -= csteps
            if self.report_every > 0 and remaining % self.report_every == 0:
                self.report(i, remaining, system)
            # check convergence every is guaranteed to check AT LEAST as often as it says
            if self.minimizer.getGlobalVariableByName("converged") == 1:
                break

        system.set_velocities(vels)
        new_state = system._context.getState(energy=True)
        new_energy = new_state.getPotentialEnergy()
        if self.max_retry > retries_so_far and new_energy > before_energy:
            with open(f"{self.sim_name}_minimization{i}.log", "a") as f:
                f.write(f"\n\nENERGY WENT UP DURING MINIMIZATION!\nfrom: {before_energy} to: {new_energy}\nRETRYING, attempts left: {self.max_retry-retries_so_far}\n\n")
                # TODO smarter scaling
            self.minimizer.global_variables["step_size"] = self.minimizer.global_variables["step_size"] * 1.03
            system.set_positions(pos)
            return self.set_reactions(reactions, system, top, i, retries_so_far + 1)
        # reporters and cleanup
        # TODO find out why velocities change significantly during minimization
        for rep in self.reporters:
            rep.post_di_minimize()
        self.remaining_eq_steps = self.eq_steps
        if self.toggle_couplings:
            system.coupling(True)
            system.reinitialize()
        if self.harmonic_constraints:
            system.constraint.constraints_to_harmonic_bonds(False)
            system.reinitialize()
        if self.remaining_eq_steps > 0:
            self.integrator.setCurrentIntegrator(2)
        else:
            self.integrator.setCurrentIntegrator(0)

    def step(self, n_steps):
        if self.remaining_eq_steps > n_steps:
            self.remaining_eq_steps -= n_steps
            self.integrator.step(n_steps * self.subdivision)
        elif self.remaining_eq_steps > 0:
            self.integrator.step(self.remaining_eq_steps * self.subdivision)
            self.finish_equilibration()
            remaining = n_steps - self.remaining_eq_steps
            self.remaining_eq_steps = 0
            if remaining > 0:
                self.integrator.step(remaining)
        else:
            self.integrator.step(n_steps)

    def finish_equilibration(self):
        self.integrator.setCurrentIntegrator(0)
        if self.reporters is not None:
            for rep in self.reporters:
                rep.post_di_equilibrate()

    def get_integrator(self):
        return self.integrator

    def getStepSize(self):
        return self.dt * mm.unit.picosecond
