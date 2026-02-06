"""
Modified version, originally from OpenMM Tools
https://github.com/choderalab/openmmtools/blob/main/openmmtools/integrators.py

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


class GradientDescentMinimizationIntegrator(mm.CustomIntegrator):

    """Simple gradient descent minimizer implemented as an integrator.

    Examples
    --------

    Create a gradient descent minimization integrator.

    >>> integrator = GradientDescentMinimizationIntegrator()

    """

    def __init__(self, initial_step_size_nm=1., smoothing_factor=0.1, random_factor=0.0):
        """
        Construct a gradient descent minimization integrator.

        Parameters
        ----------
        initial_step_size_nm
        smoothing_factor

        Notes
        -----
        An adaptive step size is used.

        """

        assert initial_step_size_nm > 0., "initial step size must be larger than 0"
        assert smoothing_factor > 0. and smoothing_factor <= 1., "smoothing factor must be between 0 and 1"
        assert random_factor >= 0., "random factor must be 0 or larger"

        timestep = 0.
        super().__init__(timestep)

        self.addGlobalVariable("step_size", initial_step_size_nm)
        self.addGlobalVariable("energy_old", 0)
        self.addGlobalVariable("energy_new", 0)
        self.addGlobalVariable("delta_energy", 0)
        self.addGlobalVariable("accept", 0)
        self.addGlobalVariable("fnorm2", 0)
        self.addPerDofVariable("x_old", 0)
        self.addPerDofVariable("est_grad", 0)
        self.addPerDofVariable("movable", 0)
        self.addGlobalVariable("eta", smoothing_factor)
        self.addGlobalVariable("random_factor", random_factor)

        # Update context state.
        self.addUpdateContextState()

        # Constrain positions.
        self.addConstrainPositions()

        # Store old energy and positions.
        self.addComputeGlobal("energy_old", "energy")
        self.addComputePerDof("x_old", "x")

        # Take step.
        self.addComputePerDof("est_grad", "(1-eta)*est_grad + eta*f + random_factor*gaussian")
        self.addComputeSum("fnorm2", "est_grad^2")
        self.addComputePerDof("x", "x+movable*step_size*est_grad/sqrt(fnorm2 + delta(fnorm2))")
        self.addConstrainPositions()

        # Ensure we only keep steps that go downhill in energy.
        self.addComputeGlobal("energy_new", "energy")
        self.addComputeGlobal("delta_energy", "energy_new-energy_old")
        # Accept also checks for NaN
        self.addComputeGlobal("accept", "step(-delta_energy) * delta(energy - energy_new)")

        self.addComputePerDof("x", "accept*x + (1-accept)*x_old")

        # Update step size.
        self.addComputeGlobal("step_size", "step_size * (2.0*accept + 0.5*(1-accept))")


