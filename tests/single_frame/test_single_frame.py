#!/usr/bin/env python3

import math
import os

import numpy as np
import openmm as mm
import openmm.app as mmapp
import pytest

from martini_daemon import Simulation, read_geometry

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
# CNAP:
# the only force in it is a 180 degree improper dihedral
# apparently there is a small relative deviation that's slightly larger than
# tolerance, but a very small absolute difference and the forces seem to pass
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
etol_override = {"CNAP": 3e-5}
ftol_override = {"cmap": 1e-2, "cutoff_LJ": 0, "pairs": 5e-5}

cutoff_nm = 1.1


@pytest.fixture
def rootdir(request):
    return os.path.dirname(request.path)


tests = [
    # notable soft skips
    "CNAP",
    "cutoff_LJ",
    # cmap, pairs
    "cmap",
    "pairs",
    "pairs_VW",
    "pairs_VWQ",
    "pairs_type",
    # biomolecule tests
    "trypsin",
    "posres",
    # polymer tests
    "polyurethane",
    # small molecule tests
    "NMC",
    "CHOL_in_W",
    "AEA",
    "NAPH",
    "nacl+waterbox",
    "solvent_mixture",
    "waterbox",
    "DPPC_DIPC_in_W",
    "CAFF",
    "BDT",
    "BZTF_CLPR",
    "benzbox",
    # specific interaction tests
    "vsite1",
    "morse",
    "vsite4fdn",
    "quartic_angle",
    "proper_dihedral",
    "cross_bond_bond",
    "vsiten2",
    "urey_bradley",
    "fourier_dihedral",
    "linear_angle",
    "fene",
    "cubic",
    "distance_restraint",
    "connection",
    "rbtorsion",
    "vsiten3",
    "vsite2fd",
    "cross_bond_angle",
    "vsite3fd",
    "g96_bond",
    "restricted_dihedral",
    "restricted_angle",
    "combined_bending_torsion",
]


# == TEST CLASS ==
class TestSingleFrame:
    def apply_constraints(self):
        """Applies constraints and vsites and checks for position change"""
        _, reference, _ = read_geometry(self.gro)
        _, respos, _ = read_geometry(self.respos)
        sim = Simulation(
            self.top, self.gro, 0, [], options={"respos": respos}, platform="Reference"
        )
        sim.context.apply_constraints()
        new_pos, box = sim.context.get_positions()
        for i in range(len(reference)):
            r_diff = box.distance(reference[i], new_pos[i])
            assert r_diff < r_tol, (
                f"Constraint/VSite position moved by {r_diff} nm (particle {i})."
            )

    def compare_daemon_gmx(self):
        _, respos, _ = read_geometry(self.respos)
        sim = Simulation(
            self.top, self.gro, 0, [], options={"respos": respos}, platform="Reference"
        )

        _, energy, _ = sim.context.get_energies()
        forces = sim.context.get_forces().flatten()
        for i in range(sim.system.num_atoms()):
            if sim.system.get_mass(i) == 0.0:
                forces[i * 3] = 0.0
                forces[i * 3 + 1] = 0.0
                forces[i * 3 + 2] = 0.0

        if energy != 0.0:
            e_diff = math.fabs(self.gmx_energy / energy - 1)
        else:
            assert self.gmx_energy == energy, f"{self.gmx_energy} != {energy}"
            e_diff = 0
        e_percent = e_diff * 100
        c_etol = etol_override.get(self.test_name) or e_tol
        assert e_diff < c_etol, (
            f"Gmx and daemon energy different by {e_percent:.2f} %.\n"
            f"Gromacs energy: {self.gmx_energy:.10e}\n"
            f"Daemon energy: {energy:.10e}\n"
            f"Relative difference {e_diff:.3e} above tolerance {c_etol:.2e}"
        )

        # Secondary comparison with raw openmm but martini daemon as parser
        # main purpose of this part is to just run the code in the public interface for exporting
        sys = sim.system.get_openmm_system()
        top = sim.get_openmm_topology()

        sim = mmapp.Simulation(
            top,
            sys,
            mm.VerletIntegrator(0.1),
            mm.Platform.getPlatformByName("Reference"),
        )
        _, pos, _ = read_geometry(self.gro)
        # note: there might be deviations here because martini_daemon's set_positions also makes vsites and constraints
        # whole across the pbc. currently no such tests exist, but in the future it may cause problems if such a test
        # is introduced.
        sim.context.setPositions(pos)
        raw_openmm_energy = (
            sim.context.getState(energy=True)
            .getPotentialEnergy()
            .value_in_unit_system(mm.unit.md_unit_system)
        )
        assert np.isclose(raw_openmm_energy, energy), "Raw OpenMM energy mismatch"

        c_ftol = ftol_override.get(self.test_name)
        if c_ftol == 0:
            return
        if c_ftol is None:
            c_ftol = f_tol

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
        box, pos, _ = read_geometry(self.gro)
        for other_atom in range(len(pos)):
            dist = box.distance(pos[atom_index], pos[other_atom])
            if np.isclose(dist, cutoff_nm):
                print(
                    f"Atoms {atom_index + 1} and {other_atom}+1 are exactly cutoff apart!"
                )
                print("This can cause artifacts in forces.")
        assert np.allclose(self.gmx_forces, forces, c_ftol, 0), (
            f"Gmx and daemon forces different by {f_percent:.2f} %.\n"
            f"Particle {atom_index + 1} (<-- indexes start from 1) "
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
        assert os.path.isfile("forces.xvg"), f"./gmxrun.sh failure for {x} (F)"
        with open("energy.xvg") as f:
            lines = f.readlines()
            self.gmx_energy = float(lines[-1].split()[-1])
        with open("forces.xvg") as f:
            lines = f.readlines()
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
        os.remove("forces.xvg")
