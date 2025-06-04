from __future__ import annotations
import openmm as mm
from .nonbonded import NonBonded


# mie math helpers
def mie_potential(x, m, n, sigma, epsilon):
    C = (n) / (n-m) * ((n/m) ** (m/(n-m)))
    return C * epsilon * ((sigma / x) ** n - (sigma / x) ** m)


def mie_potential_minimum(m, n, sigma):
    return ((m * sigma ** m) / (n * sigma ** n)) ** (1/(m-n))


def mie_force(x, m, n, sigma, epsilon):
    # Technically a misnomer, this is just V'(x),
    # a force with the right direction too would be -V'(x).
    C = (n) / (n-m) * ((n/m) ** (m/(n-m)))
    return C * epsilon * (
        -n * (sigma ** n) / (x ** (n+1)) +
        m * (sigma ** m) / (x ** (m+1))
    )


def shift(x, b, c, d, s):
    x = x - s
    return b * (x ** 2) / 2 + c * (x ** 3) / 6 + d * (x ** 4) / 24


def mie_shifted(x, m, n, b, c, d, s, sigma, epsilon):
    return mie_potential(x, m, n, sigma, epsilon) + shift(x, b, c, d, s)


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
            "step(rcut-r)*(F + S + ES);"
            f"F={C}*epsilon(type1,type2)*(term^{mie_n}-term^{mie_m});"
            "term=sigma(type1,type2)/r;"
            "S=step(switch(type1,type2)-r)*r*r*(b(type1,type2)/2+r*c(type1,type2)/6+r*r*d(type1,type2)/24);"
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

        # NOTE:
        # top_parser does sigma, epsilon -> C6, C12 parsing if non bonded type
        # is default, otherwise the values are actuall sigma and epsilon

        # add LJ parameters to the system
        C6 = []
        C12 = []
        # i,j => type index; t1,t2 => type names
        n = len(self._used_atom_types)
        for t1, i in self._used_atom_types.items():
            for t2, j in self._used_atom_types.items():
                c6, c12 = self._nb_types.get((t1, t2)) or\
                            self._nb_types.get((t2, t1)) or (0., 0.)
                C6.append(c6)
                C12.append(c12)
        # mie_n, mie_m defined above conditionally with the same condition
        # adjusted sigmas
        sigma = []
        # switch+shift param quadratic
        b = []
        # switch+shift param cubic
        c = []
        # switch+shift param quartic
        d = []
        # switch point
        switch = []
        for i in range(len(C6)):
            r_min_n_m = mie_potential_minimum(mie_m, mie_n, C6[i])
            r_min = mie_potential_minimum(6, 12, C6[i])
            sigma.append(C6[i] * r_min / r_min_n_m)

            # switching constant calculation
            F_cutoff = mie_potential(cutoff, mie_m, mie_n, sigma[i], C12[i])
            F_prime_cutoff = mie_force(cutoff, mie_m, mie_n, sigma[i], C12[i])
            diff = cutoff - r_min
            # quadratic-cubic shift+switch:
#                b.append(2 * F_prime_cutoff / diff - 6 * F_cutoff / (diff ** 2))
#                c.append(12 * F_cutoff / (diff ** 3) - 6 * F_prime_cutoff / (diff ** 2))
            # cubic-quartic shift+switch:
            b.append(0.)
            c.append(6 * F_prime_cutoff / diff ** 2 - 24 * F_cutoff / diff ** 3)
            d.append(72 / diff ** 4 * F_cutoff - 24 * F_prime_cutoff / diff ** 3)
            switch.append(r_min)
            # quadratic-cubic shift only
            # b.append(2 * F_prime_cutoff / cutoff - 6 * F_cutoff / (cutoff ** 2))
            # c.append(12 * F_cutoff / (cutoff ** 3) - 6 * F_prime_cutoff / (cutoff ** 2))
            # d.append(0.)
            # switch.append(0.)
        assert all([n * n == len(x) for x in [sigma, C12, b, c, d, switch]])
        self._force_obj.addTabulatedFunction(
            "sigma", mm.Discrete2DFunction(n, n, sigma)
        )
        self._force_obj.addTabulatedFunction(
            "epsilon", mm.Discrete2DFunction(n, n, C12)
        )
        self._force_obj.addTabulatedFunction(
            "b", mm.Discrete2DFunction(n, n, b)
        )
        self._force_obj.addTabulatedFunction(
            "c", mm.Discrete2DFunction(n, n, c)
        )
        self._force_obj.addTabulatedFunction(
            "d", mm.Discrete2DFunction(n, n, d)
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
