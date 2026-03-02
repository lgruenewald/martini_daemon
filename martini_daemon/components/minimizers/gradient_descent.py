"""
Modified version, originally from OpenMM Tools
https://github.com/choderalab/openmmtools/blob/main/openmmtools/integrators.py

Original license:
MIT License

Copyright (c) 2015-2019 Chodera lab // Memorial Sloan Kettering Cancer Center

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import openmm as mm
import numpy as np


class GradientDescentMinimizationIntegrator(mm.CustomIntegrator):

    def __init__(self, initial_step_size_nm=0.1, etol=0.0, smoothing_factor=0.1):
        """
        Construct a gradient descent minimization integrator.

        initial_step_size_nm
        Only matters at the start.
        An adaptive step size is used.

        etol
        energy tolerance, will stop doing anything once the change in energy reaches this for 1 step
        0 -> no convergence check

        smoothing_factor - 0 to 1
        the smaller the smoother but slower convergence
        """

        assert initial_step_size_nm > 0., "initial step size must be larger than 0"
        assert smoothing_factor > 0. and smoothing_factor <= 1., "smoothing factor must be between 0 and 1"

        timestep = 0.
        super().__init__(timestep)

        self.global_variables = {
            "step_size": initial_step_size_nm,
            "energy_old": 0,
            "energy_new": 0,
            "delta_energy": 0,
            "accept": 0,
            "fnorm2": 0,
            "eta": smoothing_factor,
            "probability": 0,
            "converged": 0,
            "etol": etol,
            "x_sum": 0,
            "x_sum2": 0,
            "v_sum": 0
        }

        self.per_dof_variables = {
            "x_old": 0,
            "v_old": 0,
            "est_grad": 0,
            "movable": 0
        }

        for k, v in self.global_variables.items():
            self.addGlobalVariable(k, v)

        for k, v in self.per_dof_variables.items():
            self.addPerDofVariable(k, v)

        if etol > 0:
            self.beginIfBlock("converged < 1")

        # Update context state (=Virt Sites and coupling).
        # Constrain positions.
        self.addUpdateContextState()
        self.addConstrainPositions()

        # Store old energy and positions.
        self.addComputeGlobal("energy_old", "energy")
        self.addComputePerDof("x_old", "x")
        self.addComputePerDof("v_old", "v")
        self.addComputeSum("x_sum", "x")
        self.addComputeSum("x_sum2", "x*x")
        self.addComputeSum("v_sum", "v")


        # Take step, re-constraint positions.
        self.addComputePerDof("est_grad", "(1-eta)*est_grad + eta*f")
        self.addComputeSum("fnorm2", "est_grad^2")
        self.addComputePerDof("x", "x+movable*step_size*est_grad/sqrt(fnorm2 + delta(fnorm2))")
        # x was changed, we need to do this again to actually be able to
        # see the energy
        self.addUpdateContextState()
        self.addConstrainPositions()

        # Ensure we only keep steps that go downhill in energy.
        self.addComputeGlobal("energy_new", "energy")
        self.addComputeGlobal("delta_energy", "energy_new-energy_old")

        # Accept also checks for NaN
        self.addComputeGlobal("accept", "step(-delta_energy) * delta(energy - energy_new)")

        self.beginIfBlock("accept = 0")
        # revert DoF positions
        self.addComputePerDof("x", "x_old")
        # recalc vsites to old position
#        self.addUpdateContextState()
        self.endBlock()
        self.addComputePerDof("v", "v_old")
        #self.addComputePerDof("x", "accept*x + (1-accept)*x_old")

        # Update step size.
        self.addComputeGlobal("step_size", "step_size * (2.0*accept + 0.5*(1-accept))")

        if etol > 0:
            # check convergence - must be not NaN and delta_energy < 0 and delta_energy > -etol
            self.addComputeGlobal("converged", "delta(energy-energy_new) * step(-delta_energy) * step(delta_energy + etol)")
            self.endBlock()

    def reset(self, shape):
        for k, v in self.global_variables.items():
            self.setGlobalVariableByName(k, v)

        for k, v in self.per_dof_variables.items():
            assert v == 0
            vals = np.zeros(shape)
            self.setPerDofVariableByName(k, vals)

    def report(self) -> str:
        return "\n".join(
            f"{k}: {self.getGlobalVariableByName(k)}"
            for k in self.global_variables.keys()
        ) + "\n"
