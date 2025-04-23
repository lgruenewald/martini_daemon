#!/usr/bin/env python3

import os
import openmm as mm
import openmm.app as mmapp
from openmm.unit import femtosecond, kilojoule_per_mole, kilojoule, mole, nanometer
import numpy as np
import math
import pytest

from martini_daemon.top_parser import DaemonTopFile

# == CONFIG ==
e_tol = 1e-5  # energy relative tolerance
f_tol = 1e-5  # force relative tolerance
r_tol = 2e-3  # distance tolerance


@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


tests = [
    "proper_dihedral",
    # biomolecule tests
    "trypsin",
    # polymer tests
    "polyurethane",
    # small molecule tests
    "NMC", "CHOL_in_W", "AEA", "NAPH", "nacl+waterbox", "solvent_mixture",
    "waterbox", "DPPC_DIPC_in_W", "CAFF", "BDT", "BZTF_CLPR", "benzbox",
    # specific interaction tests
    "vsite1", "pairs", "pairs_VW", "pairs_VWQ", "pairs_type",
    "morse", "vsite4fdn", "quartic_angle",
    "cross_bond_bond", "vsiten2", "urey_bradley", "fourier_dihedral",
    "linear_angle", "fene", "cubic", "distance_restraint",
    "connection", "rbtorsion", "vsiten3", "vsite2fd",
    "cross_bond_angle", "vsite3fd", "g96_bond",
    "restricted_dihedral", "restricted_angle", "combined_bending_torsion",
]


# == TEST CLASS ==
class TestSingleFrame():

    def apply_constraints(self):
        # applies constraints and vsites and checks for position change
        platform = mm.Platform.getPlatformByName("Reference")
        system, top = DaemonTopFile(self.top)
        gro = mmapp.GromacsGroFile(self.gro)
        oldpos = gro.getPositions(True)
        system.build_context(mm.VerletIntegrator(20 * femtosecond),
                             gro.getPeriodicBoxVectors(),
                             platform=platform)
        system.set_positions(gro.getPositions(True))
        system.apply_constraints()
        newpos = system.get_state().getPositions(asNumpy=True)
        r_diff = np.linalg.norm(newpos - oldpos, axis=1)
        largest_index = np.argmax(r_diff)
        nm = r_diff[largest_index]
        assert not np.any(r_diff > r_tol), (
            f"Constraint/VSite position moved by {nm} nm "
            f"(particle {largest_index})."
        )

    def compare_daemon_gmx(self):
        platform = mm.Platform.getPlatformByName("Reference")
        system, top = DaemonTopFile(self.top)
        gro = mmapp.GromacsGroFile(self.gro)
        system.build_context(mm.VerletIntegrator(20 * femtosecond),
                             gro.getPeriodicBoxVectors(),
                             platform=platform)
        system.set_positions(gro.getPositions(True))
        state = system.get_state()
        energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        forces = state.getForces(asNumpy=True).\
            value_in_unit(kilojoule / nanometer / mole).flatten()
        for vsite in system.vsites:
            forces[vsite * 3] = 0.
            forces[vsite * 3 + 1] = 0.
            forces[vsite * 3 + 2] = 0.

        e_diff = math.fabs(self.gmx_energy / energy - 1)
        e_percent = e_diff * 100
        assert e_diff < e_tol, (
            f"Gmx and daemon energy different by {e_percent:.2f} %.\n"
            f"Gromacs energy: {self.gmx_energy:.10e}\n"
            f"Daemon energy: {energy:.10e}\n"
            f"Relative difference {e_diff:.3e} above tolerance {e_tol:.2e}"
        )

        f_diff = np.fabs(self.gmx_forces - forces) / (np.fabs(forces) + f_tol)
        i_max = np.argmax(f_diff)
        max = f_diff[i_max]
        f_percent = max * 100
        force = forces[i_max]
        gmx_force = self.gmx_forces[i_max]
        abs_diff = np.fabs(force - gmx_force)
        atom_index = i_max // 3
        atom_dim = i_max % 3
        assert np.allclose(self.gmx_forces, forces, f_tol, 0), (
            f"Gmx and daemon forces different by {f_percent:.2f}.\n"
            f"Particle {atom_index} "
            f"dimension {atom_dim}\n"
            f"Absolute diff: {abs_diff:.3e}    relative diff: {max:.3e}\n"
            f"Daemon force: {force:.10e}\n"
            f"Gromacs force: {gmx_force:.10e}"
        )

    @pytest.mark.parametrize("x", tests)
    def test_single_frame(self, x, rootdir):
        os.chdir(rootdir)
        assert os.path.isfile("gmxrun.sh")
        assert os.path.isdir(x)
        os.chdir(x)
        os.system("../gmxrun.sh")
        assert os.path.isfile("energy.xvg"), f"./gmxrun.sh failure for {x} (E)"
        assert os.path.isfile("forces.xvg"), f"./gmxrun.sh failure for {x} (F)"
        with open("energy.xvg") as f:
            lines = [line for line in f]
            self.gmx_energy = float(lines[-1].split()[-1])
        with open("forces.xvg") as f:
            lines = [line for line in f]
            gmx_force_line = lines[-1].split()
            self.gmx_forces = np.array([float(x) for x in gmx_force_line][1:])
        self.top = "system.top"
        self.gro = "system.gro"
        self.apply_constraints()
        self.compare_daemon_gmx()
        os.remove("energy.xvg")
        os.remove("forces.xvg")
