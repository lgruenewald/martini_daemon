"""
Modified, originally from OpenMM Tools
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
from openmm import unit
import numpy as np


class FIREMinimizationIntegrator(mm.CustomIntegrator):
    """Fast Internal Relaxation Engine (FIRE) minimization.

    Notes
    -----
    This integrator is taken verbatim from Peter Eastman's example appearing in the CustomIntegrator header file documentation.

    References
    ----------
    Erik Bitzek, Pekka Koskinen, Franz Gaehler, Michael Moseler, and Peter Gumbsch.
    Structural Relaxation Made Simple. PRL 97:170201, 2006.
    http://dx.doi.org/10.1103/PhysRevLett.97.170201

    Examples
    --------
    >>> from openmmtools import testsystems
    >>> import openmm
    >>> t = testsystems.AlanineDipeptideVacuum()
    >>> system, positions = t.system, t.positions

    Create a FIRE integrator with default parameters and minimize for 100 steps

    >>> integrator = FIREMinimizationIntegrator()
    >>> context = openmm.Context(system, integrator)
    >>> context.setPositions(positions)
    >>> integrator.step(100)
    """

    def __init__(self, timestep_ps, tolerance=None, alpha=0.1, f_inc=1.1, f_dec=0.5, f_alpha=0.99, N_min=5):
        r"""Construct a Fast Internal Relaxation Engine (FIRE) minimization integrator.
        Parameters
        ----------
        timestep : integration timestep, in picoseconds
        tolerance : unit.Quantity compatible with kilojoules_per_mole/nanometer, optional, default = None
            Minimization will be terminated when RMS force reaches this tolerance.
        alpha : float, optional default = 0.1
            Velocity relaxation parameter, alpha \in (0,1).
        f_inc : float, optional, default = 1.1
            Timestep increment multiplicative factor.
        f_dec : float, optional, default = 0.5
            Timestep decrement multiplicative factor.
        f_alpha : float, optional, default = 0.99
            alpha multiplicative relaxation parameter
        N_min : int, optional, default = 5
            Limit on number of timesteps P is negative before decrementing timestep.
        Notes
        -----
        Velocities should be set to zero before using this integrator.
        """

        # Check input ranges.
        if not ((alpha > 0.0) and (alpha < 1.0)):
            raise Exception("alpha must be in the interval (0,1); specified alpha = %f" % alpha)

        if tolerance is None:
            tolerance = 0 * unit.kilojoules_per_mole / unit.nanometers

        timestep = timestep_ps * unit.picosecond
        dt_max = timestep * 10
        super().__init__(timestep)

        # Use high-precision constraints
        self.setConstraintTolerance(1.0e-8)

        self.global_variables = {
            "alpha": alpha,
            "P": 0,
            "N_neg": 0.,
            "fmag": 0,
            "fmax": 0,
            "ndof": 0,
            "ftol": tolerance.value_in_unit_system(unit.md_unit_system),
            "vmag": 0,
            "converged": 0,
            "E0": 0,
            "dE": 0,
            "restart": 0,
            "sanity": 0,
            "delta_t": timestep.value_in_unit_system(unit.md_unit_system)
        }

        self.per_dof_variables = {
            "x0": 0,
            "v0": 0,
            "x1": 0,
            "movable": 0
        }

        for k, v in self.global_variables.items():
            self.addGlobalVariable(k, v)

        for k, v in self.per_dof_variables.items():
            self.addPerDofVariable(k, v)

        # sanity = 0 -> finished normally
        # sanity = 1 -> in progress
        # sanity = 2 -> bad
        self.beginIfBlock("sanity = 1")
        self.addComputeGlobal("sanity", "2")
        self.endBlock()
        self.beginIfBlock("sanity = 0")
        self.addComputeGlobal("sanity", "1")
        self.endBlock()

        # Update context state.
        self.addUpdateContextState()

        # Assess convergence
        # TODO: Can we more closely match the OpenMM criterion here?
        self.beginIfBlock('converged < 1')

        # Compute fmag = |f|
        #self.addComputeGlobal('fmag', '0.0')
        self.addComputeSum('fmag', 'f*f')
        self.addComputeGlobal('fmag', 'sqrt(fmag)')

        # Compute ndof
        self.addComputeSum('ndof', '1')

        self.addComputeSum('converged', 'step(ftol - fmag/ndof)')
        self.endBlock()

        # Enclose everything in a block that checks if we have already converged.
        self.beginIfBlock('converged < 1')

        # Store old positions and energy
        self.addComputePerDof('x0', 'x')
        self.addComputePerDof('v0', 'v')
        self.addComputeGlobal('E0', 'energy')

        # MD: Take a velocity Verlet step.
        self.addComputePerDof("v", "v+0.5*delta_t*f/m")
        self.addComputePerDof("x", "x+movable*delta_t*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*delta_t*f/m+(x-x1)/delta_t")
        self.addConstrainVelocities()

        self.addComputeGlobal('dE', 'energy - E0')

        # Compute fmag = |f|
        #self.addComputeGlobal('fmag', '0.0')
        self.addComputeSum('fmag', 'f*f')
        self.addComputeGlobal('fmag', 'sqrt(fmag)')
        # Compute vmag = |v|
        #self.addComputeGlobal('vmag', '0.0')
        self.addComputeSum('vmag', 'v*v')
        self.addComputeGlobal('vmag', 'sqrt(vmag)')

        # F1: Compute P = F.v
        self.addComputeSum('P', 'f*v')

        # F2: set v = (1-alpha) v + alpha \hat{F}.|v|
        # Update velocities.
        # TODO: This must be corrected to be atomwise redirection of v magnitude along f
        self.addComputePerDof('v', '(1-alpha)*v + alpha*(f/fmag)*vmag')

        # Back up if the energy went up, protecing against NaNs
        self.addComputeGlobal('restart', '1')
        self.beginIfBlock('dE < 0')
        self.addComputeGlobal('restart', '0')
        self.endBlock()
        self.beginIfBlock('restart > 0')
        self.addComputePerDof('v', 'v0')
        self.addComputePerDof('x', 'x0')
        self.addComputeGlobal('P', '-1')
        self.endBlock()

        # If dt goes to zero, signal we've converged!
        dt_min = 1.0e-5 * timestep
        self.beginIfBlock('delta_t <= %f' % dt_min.value_in_unit_system(unit.md_unit_system))
        self.addComputeGlobal('converged', '1')
        self.endBlock()

        # F3: If P > 0 and the number of steps since P was negative > N_min,
        # Increase timestep dt = min(dt*f_inc, dt_max) and decrease alpha = alpha*f_alpha
        self.beginIfBlock('P > 0')
        # Update count of number of steps since P was negative.
        self.addComputeGlobal('N_neg', 'N_neg + 1')
        # If we have enough steps since P was negative, scale up timestep.
        self.beginIfBlock('N_neg > %d' % N_min)
        self.addComputeGlobal('delta_t', 'min(delta_t*{:f}, {:f})'.format(f_inc, dt_max.value_in_unit_system(unit.md_unit_system))) # TODO: Automatically convert dt_max to md units
        self.addComputeGlobal('alpha', 'alpha * %f' % f_alpha)
        self.endBlock()
        self.endBlock()

        # F4: If P < 0, decrease the timestep dt = dt*f_dec, freeze the system v=0,
        # and set alpha = alpha_start
        self.beginIfBlock('P < 0')
        self.addComputeGlobal('N_neg', '0.0')
        self.addComputeGlobal('delta_t', 'delta_t*%f' % f_dec)
        self.addComputePerDof('v', '0.0')
        self.addComputeGlobal('alpha', '%f' % alpha)
        self.endBlock()

        # Close block that checks for convergence.
        self.endBlock()

        # sanity back to 0
        self.beginIfBlock("sanity = 1")
        self.addComputeGlobal("sanity", "0")
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
