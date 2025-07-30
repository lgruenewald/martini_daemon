from .force import Force
import openmm as mm


class CustomReactive(Force):
    _members = 3
    _pbc = False

    def __init__(self, sysstar):
        super(CustomReactive, self).__init__(sysstar)
        # atom id -> donor id, d1 only
        self.donor_map: dict[int, int] = {}
        # atom id -> acceptor id, a1 only
        self.acceptor_map: dict[int, int] = {}
        self.reactive_types: dict[str, tuple] = {}

    def _set_force_obj(self):
        cutoff = self._sysstar.nonbonded_force.cutoff_nm
        # currently trying to model the reactive martini VSite based setup
        # but more cleanly through this

        # assumptions:
        # each SH has is a donor and acceptor at the same time -> double interactions
        # -> each d/a pair is half the energy
        # a1/d1 = SH
        # a2/d2 = COG - average of other SH and C
        # a3/d3 = C bead

        self._force_obj = mm.CustomHbondForce(
            "V_bond+V_angle+V_dihedral;"
            # dihedral
            "V_dihedral=-scaling*V_dihedral*exp(-(phi-phi0)^2/(2*w_dihedral^2));"
            "phi=abs(dihedral(a2, a1, d1, d2));"
            # angle
            "V_angle=scaling*k_angle*(theta-theta0)^2;"
            "theta=angle(a1, d1, d2);"
            # scale factor for angles and dihedrals - only should apply
            # force when bonded
            "scaling=exp(-(r-R_BB)^2/(2*w_well^2));"
            # bond
            "V_bond=barrier+step(r-R_BB)*well+step(R_BB-r)*harmonic+martini;"
            "harmonic=-V_depth+k_harmonic*(r-R_BB)^2/2;"
            "well=-V_depth*exp(-(r-R_BB)^2/(2*w_well^2));"
            "barrier=V_height * exp(-(r-R_barrier)^2/(2*w_barrier^2));"
            "martini=step(R_min-r)*(V_martini_min-V_martini);"
            "V_martini=4*epsilon*((sigma/r)^12-(sigma/r)^6);"
            "V_martini_min=4*epsilon*((sigma/R_min)^12-(sigma/R_min)^6);"
            "R_min=sigma*sqrt(2);"
            "r=distance(a1, d1);"
        )
        self._force_obj.setNonbondedMethod(
            self._force_obj.CutoffPeriodic
        )
        self._force_obj.setCutoffDistance(cutoff)
        # parameters - martini
        self._force_obj.addPerDonorParameter(
            "sigma"  # , 0.41
        )
        self._force_obj.addPerDonorParameter(
            "epsilon"  # , 2.6
        )
        # parameters - reaction barrier
        self._force_obj.addPerDonorParameter(
            "R_barrier"  # , 0.4
        )
        self._force_obj.addPerDonorParameter(
            "V_height"  # , 6
        )
        self._force_obj.addPerDonorParameter(
            "w_barrier"  # , 0.02
        )
        # parameters - reactive bond well
        self._force_obj.addPerDonorParameter(
            "R_BB"  # , 0.28
        )
        self._force_obj.addPerDonorParameter(
            "V_depth"  # , 65
        )
        self._force_obj.addPerDonorParameter(
            "w_well"  # , 0.033
        )
        # parameters - harmonic potential at short distance
        self._force_obj.addPerDonorParameter(
            "k_harmonic"  # , 250000
        )
        # parameters - angle
        self._force_obj.addPerDonorParameter(
            "theta0"  # , 127. * math.pi / 180.
        )
        self._force_obj.addPerDonorParameter(
            "k_angle"  # , 50.
        )
        # parameters - dihedral
        self._force_obj.addPerDonorParameter(
            "phi0"  # , 60. * math.pi / 180.
        )
        self._force_obj.addPerDonorParameter(
            "V_dihedral"  # , 10.
        )
        self._force_obj.addPerDonorParameter(
            "w_dihedral"  # , 0.06
        )

    def _add_to_force_obj(self, params):
        i, j, k, group_type, reactive_type = params
        ps = self.reactive_types[reactive_type]
        if group_type in {"d", "s"}:
            d = self._force_obj.addDonor(i, j, k, [*ps])
            if self.donor_map.get(i) is not None:
                raise ValueError("Only one donor per atom")
            self.donor_map[i] = d
        if group_type in {"a", "s"}:
            a = self._force_obj.addAcceptor(i, j, k, [])
            assert a == d  # always same number of acceptors and donors
            if self.acceptor_map.get(i) is not None:
                raise ValueError("Only one acceptor per atom")
            self.acceptor_map[i] = a
        if group_type == "s":
            # overlapping exclusion
            self._force_obj.addExclusion(d, a)

        # cursed exclusion logic that should be revamped later
        for pi, pj in filter(lambda x: x is not None, self._sysstar.exclusions._list):
            if pi == i:
                if self.donor_map.get(pj) is not None:
                    if group_type in {"a", "s"}:
                        self._force_obj.addExclusion(
                            self.donor_map[pj], a
                        )
                    if group_type in {"d", "s"}:
                        self._force_obj.addExclusion(
                            d, self.acceptor_map[pj]
                        )
            if pj == i:
                if self.donor_map.get(pi) is not None:
                    if group_type in {"a", "s"}:
                        self._force_obj.addExclusion(
                            self.donor_map[pi], a
                        )
                    if group_type in {"d", "s"}:
                        self._force_obj.addExclusion(
                            d, self.acceptor_map[pi]
                        )

    _filters = {"custom_reactive"}
