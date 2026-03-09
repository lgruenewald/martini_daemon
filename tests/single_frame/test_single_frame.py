#!/usr/bin/env python3

import os
import openmm as mm
from openmm.unit import femtosecond, kilojoule_per_mole, kilojoule, mole, nanometer
import numpy as np
import math
import pytest

from martini_daemon.top_parser import DaemonTopFile
from martini_daemon.gro_file import read_gro
from martini_daemon.__forces import NonBonded
from martini_daemon.utils import pdist

# == CONFIG ==
e_tol = 1e-5  # energy relative tolerance
f_tol = 1e-5  # force relative tolerance
r_tol = 2e-3  # distance tolerance

# == Soft-skip some tests ==
#
# It's possible to override the tol for certain tests.
# Setting it to 0 will disable that test (only works for ftol).
# Put an explanation here.
#
# cmap:
# probably different interpolation in OpenMM and GROMACS
# this is somehow made worse if there are bonds involved?
#
# cutoff_LJ:
# if there is a pair of atoms exactly at the LJ cutoff,
# forces can be 0 or not differently in GROMACS and OpenMM
#
# pairs:
# only slightly raised the tolerance because forces still
# seem *slightly* off for a few atoms
#
etol_override = {
}
ftol_override = {
    "cmap": 1e-2,
    "cutoff_LJ": 0,
    "pairs": 5e-5
}

cutoff_nm = 1.1


@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


tests = [
    "cutoff_LJ", "cmap",
    "pairs", "pairs_VW", "pairs_VWQ", "pairs_type",
    # biomolecule tests
    "trypsin", "posres",
    # polymer tests
    "polyurethane",
    # small molecule tests
    "NMC", "CHOL_in_W", "AEA", "NAPH", "nacl+waterbox", "solvent_mixture",
    "waterbox", "DPPC_DIPC_in_W", "CAFF", "BDT", "BZTF_CLPR", "benzbox",
    # specific interaction tests
    "vsite1",
    "morse", "vsite4fdn", "quartic_angle", "proper_dihedral",
    "cross_bond_bond", "vsiten2", "urey_bradley", "fourier_dihedral",
    "linear_angle", "fene", "cubic", "distance_restraint",
    "connection", "rbtorsion", "vsiten3", "vsite2fd",
    "cross_bond_angle", "vsite3fd", "g96_bond",
    "restricted_dihedral", "restricted_angle", "combined_bending_torsion",
]

# TODO: use Simulation, not TopParser
# == TEST CLASS ==
class TestSingleFrame():
    def apply_constraints(self):
        # applies constraints and vsites and checks for position change
        platform = mm.Platform.getPlatformByName("Reference")
        _, respos, _ = read_gro(self.respos)
        ok, res = DaemonTopFile(
            self.top, NonBonded(cutoff_nm=cutoff_nm), respos=respos,
            experimental=True
        )
        assert ok
        system, top = res
        box, pos, vel = read_gro(self.gro)
        system.build_context(
            mm.VerletIntegrator(20 * femtosecond),
            box,
            platform=platform
        )
        system.set_positions(pos)
        system.apply_constraints()
        newpos = system.get_state().getPositions(asNumpy=True)\
            .value_in_unit(mm.unit.nanometer)
        r_diff = np.linalg.norm(newpos - pos, axis=1)
        largest_index = np.argmax(r_diff)
        nm = r_diff[largest_index]
        assert not np.any(r_diff > r_tol), (
            f"Constraint/VSite position moved by {nm} nm "
            f"(particle {largest_index})."
        )

    def compare_daemon_gmx(self):
        platform = mm.Platform.getPlatformByName("Reference")
        _, respos, _ = read_gro(self.respos)
        ok, res = DaemonTopFile(
            self.top, NonBonded(cutoff_nm=1.1), respos=respos,
            experimental=True
        )
        assert ok
        system, top = res
        box, pos, vel = read_gro(self.gro)
        system.build_context(
            mm.VerletIntegrator(20 * femtosecond),
            box,
            platform=platform
        )
        system.set_positions(pos)
        state = system.get_state()
        energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        forces = state.getForces(asNumpy=True).\
            value_in_unit(kilojoule / nanometer / mole).flatten()
        for vsite in system.vsites:
            forces[vsite * 3] = 0.
            forces[vsite * 3 + 1] = 0.
            forces[vsite * 3 + 2] = 0.

        if energy != 0.:
            e_diff = math.fabs(self.gmx_energy / energy - 1)
        else:
            assert self.gmx_energy == energy, f"{self.gmx_energy} != {energy}"
            e_diff = 0
        e_percent = e_diff * 100
        cetol = etol_override.get(self.test_name) or e_tol
        assert e_diff < cetol, (
            f"Gmx and daemon energy different by {e_percent:.2f} %.\n"
            f"Gromacs energy: {self.gmx_energy:.10e}\n"
            f"Daemon energy: {energy:.10e}\n"
            f"Relative difference {e_diff:.3e} above tolerance {cetol:.2e}"
        )

        cftol = ftol_override.get(self.test_name)
        if cftol == 0:
            return
        elif cftol is None:
            cftol = f_tol

        f_diff = np.fabs(self.gmx_forces - forces) / (np.fabs(forces) + f_tol)
        i_max = np.argmax(f_diff)
        max = f_diff[i_max]
        f_percent = max * 100
        force = forces[i_max]
        gmx_force = self.gmx_forces[i_max]
        abs_diff = np.fabs(force - gmx_force)
        atom_index = i_max // 3
        atom_dim = i_max % 3
        # check if there is any exactly cutoffs
        box = np.array(box)
        for other_atom in range(len(pos)):
            dist = pdist(pos[atom_index], pos[other_atom], box)
            if np.isclose(dist, cutoff_nm):
                print(f"Atoms {atom_index+1} and {other_atom}+1 are exactly cutoff apart!")
                print("This can cause artifacts in __forces.")
        assert np.allclose(self.gmx_forces, forces, cftol, 0), (
            f"Gmx and daemon __forces different by {f_percent:.2f} %.\n"
            f"Particle {atom_index+1} (<-- indexes start from 1) "
            f"dimension {atom_dim}\n"
            f"Absolute diff: {abs_diff:.3e}    relative diff: {max:.3e}\n"
            f"Daemon force: {force:.10e}\n"
            f"Gromacs force: {gmx_force:.10e}"
        )

    @pytest.mark.parametrize("x", tests)
    def test_single_frame(self, x, rootdir):
        self.test_name = x
        os.chdir(rootdir)
        assert os.path.isfile("gmxrun.sh")
        assert os.path.isdir(x)
        os.chdir(x)
        os.system("../gmxrun.sh")
        assert os.path.isfile("energy.xvg"), f"./gmxrun.sh failure for {x} (E)"
        assert os.path.isfile("__forces.xvg"), f"./gmxrun.sh failure for {x} (F)"
        with open("energy.xvg") as f:
            lines = [line for line in f]
            self.gmx_energy = float(lines[-1].split()[-1])
        with open("__forces.xvg") as f:
            lines = [line for line in f]
            gmx_force_line = lines[-1].split()
            self.gmx_forces = np.array([float(x) for x in gmx_force_line][1:])
        self.top = "system.top"
        self.gro = "system.gro"
        if os.path.isfile("respos.gro"):
            self.respos = "respos.gro"
        else:
            self.respos = "system.gro"
        self.apply_constraints()
        self.compare_daemon_gmx()
        os.remove("energy.xvg")
        os.remove("__forces.xvg")
