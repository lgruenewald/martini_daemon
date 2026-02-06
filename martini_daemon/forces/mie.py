from __future__ import annotations
import openmm as mm
from .nonbonded import NonBonded
import numpy as np


def C(n, m):
    """
    C factor for Mie n, m
    """
    return (n) / (n-m) * ((n/m) ** (m/(n-m)))

def mie_potential(x, n, m, sigma, epsilon):
    """
    V(x) mie potential for orders m>n with sigma and epsilon as params.
    """
    return C(n, m) * epsilon * ((sigma / x) ** n - (sigma / x) ** m)


def mie_potential_minimum(n, m, sigma):
    """
    Returns r at which the potential is at its minimum for
    a particular n, m and sigma.
    """
    return ((m * sigma ** m) / (n * sigma ** n)) ** (1/(m-n))


def mie_derivative(x, n, m, sigma, epsilon):
    """
    V'(x) derivative of the mie potential.
    """
    return C(n, m) * epsilon * (
        -n * (sigma ** n) / (x ** (n+1))
        + m * (sigma ** m) / (x ** (m+1))
    )

def mie_2derivative(x, n, m, sigma, epsilon):
    """
    V''(x) second derivative of the mie potential.
    """
    return C(n, m) * epsilon * (
        -n * (n + 1) * (sigma ** n) / (x ** (n + 2))
        + m * (m + 1) * (sigma ** m) / (x ** (m + 2))
    )

def derivative_shifted(x, n, m, a, b, c, s, sigma, epsilon):
    """
    V'(x) shifted mie potential derivative for order m>n.
    Args:
    b, c, s - shift args
    n, m, sigma, epsilon - mie_derivative() args
    """
    return mie_derivative(x, n, m, sigma, epsilon) + (
        b * (x - s) ** 2 / 2
        + c * (x - s) ** 3 / 6
    )

def potential_shifted(x, n, m, a, b, c, s, sigma, epsilon):
    """
    V(x) shifted mie potential for order m>n.
    Args:
    b, c, s - shift() args
    n, m, sigma, epsilon - mie_potential() args
    """
    return mie_potential(x, n, m, sigma, epsilon) + (
        a
        + b * (x - s) ** 3 / 6
        + c * (x - s) ** 4 / 24
    )

def convert_parameters(n, m, sigma, epsilon, cutoff):
    """
    Given a LJ Martini sigma and epsilon and a mie order of m>n,
    returns a tuple of:
    - sigma - the new Mie sigma
    - epsilon - the new Mie epsilon
    - b, c - correction Taylor expansion constants
    - switch - switching distance
    """
    # new sigma calculation - keep minimum where it was
    r_min_mie = mie_potential_minimum(n, m, sigma)
    r_min_LJ = mie_potential_minimum(12, 6, sigma)
    new_sigma = sigma * r_min_LJ / r_min_mie

    # switching constant calculation
    Fc = mie_derivative(cutoff, n, m, new_sigma, epsilon)
    Fpc = mie_2derivative(cutoff, n, m, new_sigma, epsilon)
    diff = cutoff - r_min_LJ
    b = 2 * Fpc / diff - 6 * Fc / (diff ** 2)
    c = 12 * Fc / (diff ** 3) - 6 * Fpc / (diff ** 2)
    a = -potential_shifted(
        cutoff, n, m, 0, b, c, r_min_LJ,
        new_sigma, epsilon
    )

    # epsilon adjustment - BROKEN?
    # mie integral
    res = 1000
    x = np.linspace(new_sigma, cutoff, res)
    V = potential_shifted(x, n, m, a, b, c, r_min_LJ, new_sigma, epsilon)
    # it's not so easy writing an analytical integral for any n/m pair
    # (with switch/shift)
    integral = np.trapz(V * (x ** 2), x)
    # LJ integral
    integral_LJ = -8/9 * sigma ** 3 * epsilon
    new_epsilon = epsilon * (integral_LJ / integral)
    assert integral_LJ / integral > 0.

    return (
        new_sigma,
        new_epsilon,
        a, b, c,
        r_min_LJ
    )
    

class MiePotential(NonBonded):

    def __init__(self, mie_n, mie_m, epsilon_r=15.0, cutoff_nm=1.1):
        super(MiePotential, self).__init__(epsilon_r, cutoff_nm)
        self.mie_n = mie_n
        self.mie_m = mie_m

    def _build(self):
        mie_n = self.mie_n
        mie_m = self.mie_m
        cutoff = self.cutoff_nm
        epsilon_r = self.epsilon_r
        C = (mie_n) / (mie_n-mie_m) * ((mie_n/mie_m) ** (mie_m/(mie_n-mie_m)))
        self._force_obj = mm.CustomNonbondedForce(
            "step(rcut-r)*(F + S + a(type1,type2) + ES);"
            f"F={C}*epsilon(type1,type2)*(term^{mie_n}-term^{mie_m});"
            "term=sigma(type1,type2)/r;"
            "S=step(diff)*(b(type1,type2)/6*diff^3+c(type1,type2)/24*diff^4);"
            "diff=r-switch(type1,type2);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            "f = 138.935458;"
            f"rcut={cutoff};"
        )

        self._force_obj.addPerParticleParameter("type")
        self._force_obj.addPerParticleParameter("q")
        self._force_obj.setNonbondedMethod(
            mm.CustomNonbondedForce.CutoffPeriodic
        )
        self._force_obj.setCutoffDistance(self.cutoff_nm)

        for i in range(self._sysstar.len_atoms()):
            type, charge, _ = self._sysstar.get_atom_details(i)
            atom_type_id = self.use_atom_type(type)
            self._force_obj.addParticle([atom_type_id, charge])

        for (i, j) in filter(None, self._exclusions._list):
            self._force_obj.addExclusion(i, j)

        # add LJ parameters to the system
        lj_sigma = []
        lj_epsilon = []
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                sigma, epsilon = self._nb_types.get((t1, t2)) or\
                            self._nb_types.get((t2, t1)) or (0., 0.)
                lj_sigma.append(sigma)
                lj_epsilon.append(epsilon)
        # adjusted sigmas
        sigma = []
        # adjusted epsilons
        epsilon = []
        # shift linear
        a = []
        # switch+shift param quadratic (in force)
        b = []
        # switch+shift param cubic (in force)
        c = []
        # switch point
        switch = []
        for i in range(len(lj_sigma)):
            (
                new_sigma, new_epsilon,
                new_a, new_b, new_c,
                new_switch
            ) = convert_parameters(mie_m, mie_n, lj_sigma[i], lj_epsilon[i], cutoff)
            sigma.append(new_sigma)
            epsilon.append(new_epsilon)
            a.append(new_a)
            b.append(new_b)
            c.append(new_c)
            switch.append(new_switch)

        assert all([n * n == len(x) for x in [sigma, epsilon, b, c, switch]])
        self._force_obj.addTabulatedFunction(
            "sigma", mm.Discrete2DFunction(n, n, sigma)
        )
        self._force_obj.addTabulatedFunction(
            "epsilon", mm.Discrete2DFunction(n, n, epsilon)
        )
        self._force_obj.addTabulatedFunction(
            "a", mm.Discrete2DFunction(n, n, a)
        )
        self._force_obj.addTabulatedFunction(
            "b", mm.Discrete2DFunction(n, n, b)
        )
        self._force_obj.addTabulatedFunction(
            "c", mm.Discrete2DFunction(n, n, c)
        )
        self._force_obj.addTabulatedFunction(
            "switch", mm.Discrete2DFunction(n, n, switch)
        )

    def _additional_save(self, f):
        super(MiePotential, self)._additional_save(f)
        f.dump(self.mie_m)
        f.dump(self.mie_n)

    def _additional_load(self, f):
        super(MiePotential, self)._additional_load(f)
        self.mie_m = f.load()
        self.mie_n = f.load()
